#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stratified Epidemiological Analysis – Response to Reviewer #3 (Point 3)
========================================================================

Generates the additional stratified tables and figures requested by Reviewer #3:
  1. Mean annual incidence rate by age group and sex (2011–2023)
  2. Prevalence by age group and sex (2023)
  3. Region-specific crude and age-standardized prevalence (2023) – formatted table

Outputs
-------
  output/tables/table_incidence_by_age_sex.csv
  output/tables/table_prevalence_by_age_sex_2023.csv
  output/tables/table_regional_asr_prevalence_2023.csv
  output/figures/fig_incidence_by_age_sex.png
  output/figures/fig_prevalence_by_age_sex_2023.png
"""

import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2

from analisis_em.config import (
    EXCEL_FILE, INE_POP_FILE, PREV_REGION_FILE, BEN_SYS_FILE,
    MIN_YEAR, MAX_YEAR, INC_START_YEAR, TOP_OPEN,
    get_output_dirs,
)
from analisis_em.io_data import read_bdem, load_ine_population
from analisis_em.incidence import compute_incidence_from_bdem
from analisis_em.prevalence import build_prevalence_base
from analisis_em.utils import to_age_group_5y, age_group_sortkey, sort_age_index
from analisis_em.who_standards import who_world_standard_weights_collapsed
from analisis_em.statistics import poisson_rate_ci, asr_and_ci_from_counts

# Suppress non-critical warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)


# =============================================================================
# HELPERS
# =============================================================================

AGE_ORDER = [
    "0-4", "5-9", "10-14", "15-19", "20-24", "25-29",
    "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
    "60-64", "65-69", "70-74", "75-79", f"{TOP_OPEN}+",
]


def _sort_age(series: pd.Index) -> list:
    return sorted(series, key=age_group_sortkey)


def poisson_ci_rate(cases: int, population: float, mult: float = 1e5) -> tuple:
    """Return (rate, lcl, ucl) per 100 000 using exact Poisson CI."""
    if population <= 0 or np.isnan(population):
        return (np.nan, np.nan, np.nan)
    r, l, u = poisson_rate_ci(int(cases), float(population), mult)
    return r, l, u


def format_rate_ci(r, l, u, decimals: int = 2) -> str:
    if any(np.isnan(v) for v in [r, l, u]):
        return "–"
    return f"{r:.{decimals}f} ({l:.{decimals}f}–{u:.{decimals}f})"


# =============================================================================
# 1. MEAN ANNUAL INCIDENCE BY AGE GROUP AND SEX (2011-2023)
# =============================================================================

def table_incidence_by_age_sex(
    incident: pd.DataFrame,
    pop_long: pd.DataFrame,
) -> pd.DataFrame:
    """
    Pool all incident cases and population over 2011-2023 to compute a
    single mean-annual age- and sex-specific incidence rate per 100 000.

    Pooled rate = sum(cases_y) / sum(pop_y)  ×  100 000   for each (age, sex) cell.
    95% CI via exact Poisson on the pooled count / pooled person-years.
    """
    years = list(range(INC_START_YEAR + 1, MAX_YEAR + 1))
    inc = incident[incident["year"].isin(years)].copy()

    pop = pop_long[pop_long["year"].isin(years)].copy()
    pop["edad_grupo"] = pop["Edad"].apply(lambda x: to_age_group_5y(x, TOP_OPEN))

    # Pool population (person-years)
    pop_pool = (
        pop.groupby(["edad_grupo", "sexo"])["poblacion"].sum()
    )
    pop_pool_total = pop.groupby("edad_grupo")["poblacion"].sum()

    # Pool cases
    cases_both = inc.groupby("edad_grupo").size().rename("cases")
    cases_sex  = inc.groupby(["edad_grupo", "sexo"]).size().rename("cases")

    # All age groups observed
    age_groups = _sort_age(
        set(list(cases_both.index) + [g for g, _ in cases_sex.index])
    )

    rows = []
    for ag in age_groups:
        n_total  = int(cases_both.get(ag, 0))
        n_f      = int(cases_sex.get((ag, "F"), 0))
        n_m      = int(cases_sex.get((ag, "M"), 0))

        pop_t  = float(pop_pool_total.get(ag, np.nan))
        pop_f  = float(pop_pool.get((ag, "F"), np.nan))
        pop_m  = float(pop_pool.get((ag, "M"), np.nan))

        r_t, l_t, u_t = poisson_ci_rate(n_total, pop_t)
        r_f, l_f, u_f = poisson_ci_rate(n_f, pop_f)
        r_m, l_m, u_m = poisson_ci_rate(n_m, pop_m)

        rows.append({
            "Age group":          ag,
            "Total cases (N)":    n_total,
            "Female cases (N)":   n_f,
            "Male cases (N)":     n_m,
            "Total rate (95% CI)":   format_rate_ci(r_t, l_t, u_t),
            "Female rate (95% CI)":  format_rate_ci(r_f, l_f, u_f),
            "Male rate (95% CI)":    format_rate_ci(r_m, l_m, u_m),
            # raw numbers for plotting
            "_r_t": r_t, "_r_f": r_f, "_r_m": r_m,
            "_l_t": l_t, "_l_f": l_f, "_l_m": l_m,
            "_u_t": u_t, "_u_f": u_f, "_u_m": u_m,
        })

    return pd.DataFrame(rows)


# =============================================================================
# 2. PREVALENCE BY AGE GROUP AND SEX – 2023
# =============================================================================

def table_prevalence_by_age_sex(
    base_prev: pd.DataFrame,
    pop_long: pd.DataFrame,
    year: int = MAX_YEAR,
) -> pd.DataFrame:
    """
    Age- and sex-specific prevalence rates per 100 000 for a given year.
    """
    alive = base_prev[base_prev["death_is_null"] & (base_prev["diag_year"] <= year)].copy()

    pop_y = pop_long[pop_long["year"] == year].copy()
    pop_y["edad_grupo"] = pop_y["Edad"].apply(lambda x: to_age_group_5y(x, TOP_OPEN))

    pop_t  = pop_y.groupby("edad_grupo")["poblacion"].sum()
    pop_sx = pop_y.groupby(["edad_grupo", "sexo"])["poblacion"].sum()

    cases_t  = alive.groupby("edad_grupo").size().rename("cases")
    cases_sx = alive.groupby(["edad_grupo", "sexo"]).size().rename("cases")

    age_groups = _sort_age(
        set(list(cases_t.index) + [g for g, _ in cases_sx.index])
    )

    rows = []
    for ag in age_groups:
        n_total = int(cases_t.get(ag, 0))
        n_f     = int(cases_sx.get((ag, "F"), 0))
        n_m     = int(cases_sx.get((ag, "M"), 0))

        p_t = float(pop_t.get(ag, np.nan))
        p_f = float(pop_sx.get((ag, "F"), np.nan))
        p_m = float(pop_sx.get((ag, "M"), np.nan))

        r_t, l_t, u_t = poisson_ci_rate(n_total, p_t)
        r_f, l_f, u_f = poisson_ci_rate(n_f, p_f)
        r_m, l_m, u_m = poisson_ci_rate(n_m, p_m)

        rows.append({
            "Age group":          ag,
            "Total cases (N)":    n_total,
            "Female cases (N)":   n_f,
            "Male cases (N)":     n_m,
            "Total rate (95% CI)":   format_rate_ci(r_t, l_t, u_t),
            "Female rate (95% CI)":  format_rate_ci(r_f, l_f, u_f),
            "Male rate (95% CI)":    format_rate_ci(r_m, l_m, u_m),
            "_r_t": r_t, "_r_f": r_f, "_r_m": r_m,
            "_l_t": l_t, "_l_f": l_f, "_l_m": l_m,
            "_u_t": u_t, "_u_f": u_f, "_u_m": u_m,
        })

    return pd.DataFrame(rows)


# =============================================================================
# 3. REGIONAL ASR PREVALENCE TABLE – 2023
# =============================================================================

def table_regional_asr(outdir: Path) -> pd.DataFrame:
    """
    Read the existing prevalence_latitude_merge_2023.csv and return a
    publication-ready table with crude and ASR prevalence by region.
    """
    src = outdir / "prevalence_latitude_merge_2023.csv"
    if not src.exists():
        raise FileNotFoundError(f"Missing: {src}")

    df = pd.read_csv(src)
    df = df.rename(columns={
        "region_name":      "Region",
        "crude_prev_per_100k": "Crude prevalence (per 100,000)",
        "who_asr_per_100k": "Age-standardized prevalence (per 100,000)",
        "latitud":          "Latitude (°S)",
    })

    df["Crude prevalence (per 100,000)"] = df["Crude prevalence (per 100,000)"].round(2)
    df["Age-standardized prevalence (per 100,000)"] = df["Age-standardized prevalence (per 100,000)"].round(2)

    return df[["Region", "Crude prevalence (per 100,000)",
               "Age-standardized prevalence (per 100,000)", "Latitude (°S)"]].copy()


# =============================================================================
# FIGURES
# =============================================================================

def fig_incidence_age_sex(df: pd.DataFrame, path: Path) -> None:
    """
    Grouped bar chart: mean annual incidence rate per 100 000 by age group and sex.
    """
    df_plot = df.dropna(subset=["_r_t"]).copy()
    ag = df_plot["Age group"].tolist()
    x  = np.arange(len(ag))
    w  = 0.28

    fig, ax = plt.subplots(figsize=(13, 6))

    b_t = ax.bar(x - w, df_plot["_r_t"], w,
                 yerr=[df_plot["_r_t"] - df_plot["_l_t"],
                       df_plot["_u_t"] - df_plot["_r_t"]],
                 capsize=3, label="Total", color="#3b82f6", alpha=0.85)
    b_f = ax.bar(x,     df_plot["_r_f"], w,
                 yerr=[df_plot["_r_f"] - df_plot["_l_f"],
                       df_plot["_u_f"] - df_plot["_r_f"]],
                 capsize=3, label="Female", color="#ef4444", alpha=0.85)
    b_m = ax.bar(x + w, df_plot["_r_m"], w,
                 yerr=[df_plot["_r_m"] - df_plot["_l_m"],
                       df_plot["_u_m"] - df_plot["_r_m"]],
                 capsize=3, label="Male", color="#10b981", alpha=0.85)

    ax.set_xticks(x)
    ax.set_xticklabels(ag, rotation=45, ha="right", fontsize=10)
    ax.set_ylabel("Mean annual incidence rate per 100,000", fontsize=12)
    ax.set_xlabel("Age group (years)", fontsize=12)
    ax.set_title(
        f"Mean Annual Incidence Rate of RRMS by Age Group and Sex\n"
        f"Chile, {INC_START_YEAR + 1}–{MAX_YEAR}",
        fontsize=13,
    )
    ax.legend(fontsize=11)
    ax.grid(True, axis="y", alpha=0.35, linestyle="--")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


def fig_prevalence_age_sex(df: pd.DataFrame, path: Path, year: int = MAX_YEAR) -> None:
    """
    Grouped bar chart: prevalence rate per 100 000 by age group and sex.
    """
    df_plot = df.dropna(subset=["_r_t"]).copy()
    ag = df_plot["Age group"].tolist()
    x  = np.arange(len(ag))
    w  = 0.28

    fig, ax = plt.subplots(figsize=(13, 6))

    ax.bar(x - w, df_plot["_r_t"], w,
           yerr=[df_plot["_r_t"] - df_plot["_l_t"],
                 df_plot["_u_t"] - df_plot["_r_t"]],
           capsize=3, label="Total", color="#3b82f6", alpha=0.85)
    ax.bar(x,     df_plot["_r_f"], w,
           yerr=[df_plot["_r_f"] - df_plot["_l_f"],
                 df_plot["_u_f"] - df_plot["_r_f"]],
           capsize=3, label="Female", color="#ef4444", alpha=0.85)
    ax.bar(x + w, df_plot["_r_m"], w,
           yerr=[df_plot["_r_m"] - df_plot["_l_m"],
                 df_plot["_u_m"] - df_plot["_r_m"]],
           capsize=3, label="Male", color="#10b981", alpha=0.85)

    ax.set_xticks(x)
    ax.set_xticklabels(ag, rotation=45, ha="right", fontsize=10)
    ax.set_ylabel("Prevalence rate per 100,000", fontsize=12)
    ax.set_xlabel("Age group (years)", fontsize=12)
    ax.set_title(
        f"Prevalence of RRMS by Age Group and Sex – Chile, {year}",
        fontsize=13,
    )
    ax.legend(fontsize=11)
    ax.grid(True, axis="y", alpha=0.35, linestyle="--")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    outdir, gdir = get_output_dirs()

    print("=" * 60)
    print("LOADING DATA...")
    print("=" * 60)
    df       = read_bdem(EXCEL_FILE)
    pop_long = load_ine_population(INE_POP_FILE)

    print("=" * 60)
    print("COMPUTING INCIDENT CASES...")
    print("=" * 60)
    incident  = compute_incidence_from_bdem(df)
    base_prev = build_prevalence_base(df)

    # ------------------------------------------------------------------
    # TABLE 1: Mean annual incidence by age group and sex
    # ------------------------------------------------------------------
    print("\n--- Table: Mean annual incidence by age group & sex (2011–2023) ---")
    t_inc = table_incidence_by_age_sex(incident, pop_long)

    cols_export = [c for c in t_inc.columns if not c.startswith("_")]
    t_inc[cols_export].to_csv(
        outdir / "table_incidence_by_age_sex.csv",
        index=False, encoding="utf-8-sig"
    )
    print(t_inc[cols_export].to_string(index=False))

    # ------------------------------------------------------------------
    # TABLE 2: Prevalence by age group and sex – 2023
    # ------------------------------------------------------------------
    print(f"\n--- Table: Prevalence by age group & sex ({MAX_YEAR}) ---")
    t_prev = table_prevalence_by_age_sex(base_prev, pop_long, year=MAX_YEAR)

    t_prev[cols_export].to_csv(
        outdir / "table_prevalence_by_age_sex_2023.csv",
        index=False, encoding="utf-8-sig"
    )
    print(t_prev[cols_export].to_string(index=False))

    # ------------------------------------------------------------------
    # TABLE 3: Regional ASR (formatted from existing CSV)
    # ------------------------------------------------------------------
    print("\n--- Table: Regional ASR prevalence (2023) ---")
    t_reg = table_regional_asr(outdir)
    t_reg.to_csv(
        outdir / "table_regional_asr_prevalence_2023.csv",
        index=False, encoding="utf-8-sig"
    )
    print(t_reg.to_string(index=False))

    # ------------------------------------------------------------------
    # FIGURES
    # ------------------------------------------------------------------
    print("\n--- Generating figures ---")
    fig_incidence_age_sex(
        t_inc,
        gdir / "fig_incidence_by_age_sex.png"
    )
    fig_prevalence_age_sex(
        t_prev,
        gdir / "fig_prevalence_by_age_sex_2023.png"
    )

    print("\n" + "=" * 60)
    print("STRATIFIED ANALYSIS COMPLETE")
    print(f"Tables  -> {outdir}")
    print(f"Figures -> {gdir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
