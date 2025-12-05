"""
Visualizaciones
===============

Funciones para generación de gráficos:
- Gráficos de líneas y barras básicos
- Series temporales con bandas de confianza
- Histogramas por edad
- Gráficos comparativos
"""

from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

from .config import INC_START_YEAR, MAX_YEAR, MIN_YEAR, ISAPRE_NAMES, TOP_OPEN
from .utils import age_group_sortkey, to_age_group_5y
from .statistics import pearson_ci, ols_confidence_band, poisson_rate_ci
from .who_standards import who_world_standard_weights_collapsed


# =============================================================================
# ESTILOS COMUNES
# =============================================================================

def _poster_axes(ax):
    """Aplica estilo común a los ejes."""
    ax.grid(True, axis="y", alpha=0.35, linestyle="--", linewidth=0.8)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="both", labelsize=11)
    ax.set_ylim(bottom=0)


def _placeholder(ax, title: str, xlabel: str, ylabel: str):
    """Crea gráfico placeholder cuando no hay datos."""
    ax.set_title(title, fontsize=14)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    _poster_axes(ax)
    ax.text(0.5, 0.5, "No data available", ha="center", va="center",
            transform=ax.transAxes, fontsize=13, alpha=0.8)


# =============================================================================
# GRÁFICOS DE LÍNEAS BÁSICOS
# =============================================================================

def save_line(series: pd.Series, title: str, path: Path,
              ylabel: str = "Rate per 100,000", xlab: str = "Year",
              x_index: range | None = None) -> None:
    """
    Guarda gráfico de línea simple.
    
    Args:
        series: Datos a graficar
        title: Título del gráfico
        path: Ruta de salida
        ylabel: Etiqueta eje Y
        xlab: Etiqueta eje X
        x_index: Índice opcional para reindexar
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    s = pd.Series(series) if series is not None else pd.Series(dtype=float)
    
    if x_index is not None:
        s = s.reindex(list(x_index))
    
    if s.dropna().empty:
        _placeholder(ax, title, xlab, ylabel)
    else:
        s = s.sort_index()
        ax.plot(s.index, s.values, marker="o", linewidth=2.8)
        ax.set_title(title, fontsize=14)
        ax.set_xlabel(xlab, fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        _poster_axes(ax)
    
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def save_line_multi(series_dict: dict[str, pd.Series], title: str, path: Path,
                    ylabel: str = "Rate per 100,000", xlab: str = "Year",
                    x_index: range | None = None, legend_title: str | None = None) -> None:
    """
    Guarda gráfico de múltiples líneas.
    
    Args:
        series_dict: Dict {label: Serie}
        title: Título
        path: Ruta de salida
        ylabel, xlab: Etiquetas de ejes
        x_index: Índice opcional
        legend_title: Título de leyenda
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    has_any = False
    
    for label, ser in (series_dict or {}).items():
        s = pd.Series(ser) if ser is not None else pd.Series(dtype=float)
        if x_index is not None:
            s = s.reindex(list(x_index))
        if not s.dropna().empty:
            has_any = True
            ax.plot(s.index, s.values, marker="o", linewidth=2.4, label=label)
    
    if not has_any:
        _placeholder(ax, title, xlab, ylabel)
    else:
        ax.set_title(title, fontsize=14)
        ax.set_xlabel(xlab, fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        _poster_axes(ax)
        ax.legend(frameon=False, title=legend_title)
    
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def save_two_series_line(a: pd.Series, b: pd.Series, labels: tuple[str, str],
                         title: str, path: Path, x_index: range) -> None:
    """Guarda gráfico comparativo de dos series."""
    fig, ax = plt.subplots(figsize=(10, 6))
    sa = pd.Series(a).reindex(list(x_index)) if a is not None else pd.Series(dtype=float)
    sb = pd.Series(b).reindex(list(x_index)) if b is not None else pd.Series(dtype=float)
    
    if sa.dropna().empty and sb.dropna().empty:
        _placeholder(ax, title, "Year", "Rate per 100,000")
    else:
        if not sa.dropna().empty:
            ax.plot(sa.index, sa.values, marker="o", linewidth=2.8, label=labels[0])
        if not sb.dropna().empty:
            ax.plot(sb.index, sb.values, marker="o", linewidth=2.8, label=labels[1])
        ax.set_title(title, fontsize=14)
        ax.set_xlabel("Year")
        ax.set_ylabel("Rate per 100,000")
        _poster_axes(ax)
        ax.legend(frameon=False)
    
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


# =============================================================================
# GRÁFICOS DE BARRAS
# =============================================================================

def save_bar(series: pd.Series, title: str, path: Path,
             ylabel: str = "Count", xlab: str = "Year",
             x_index: range | None = None) -> None:
    """Guarda gráfico de barras."""
    fig, ax = plt.subplots(figsize=(10, 6))
    s = pd.Series(series) if series is not None else pd.Series(dtype=float)
    
    if x_index is not None:
        s = s.reindex(list(x_index))
    
    if s.dropna().empty:
        _placeholder(ax, title, xlab, ylabel)
    else:
        s = s.sort_index()
        s.plot(kind="bar", ax=ax)
        ax.set_title(title, fontsize=14)
        ax.set_xlabel(xlab, fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        _poster_axes(ax)
    
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


# =============================================================================
# SERIES TEMPORALES CON BANDAS DE CONFIANZA
# =============================================================================

def _plot_timeseries_with_ci(ax, frames: dict, title: str, ylabel: str = "Rate per 100,000"):
    """
    Dibuja series temporales con ribbons de IC.
    
    Args:
        ax: Ejes de matplotlib
        frames: Dict {label: DataFrame con rate, lcl, ucl}
        title: Título
        ylabel: Etiqueta eje Y
    """
    has_any = False
    
    for lbl, df in frames.items():
        if df is None or df.empty:
            continue
        df = df.sort_index()
        has_any = True
        ax.plot(df.index, df["rate"], label=lbl, linewidth=2.2, marker="o")
        ax.fill_between(df.index, df["lcl"], df["ucl"], alpha=0.18, linewidth=0)
    
    if not has_any:
        _placeholder(ax, title, "Year", ylabel)
        return
    
    ax.set_title(title, fontsize=14)
    ax.set_xlabel("Year")
    ax.set_ylabel(ylabel)
    _poster_axes(ax)
    ax.legend(frameon=False, loc="best")


def plot_incidence_by_year_crude_and_asr_with_ci(incident: pd.DataFrame, 
                                                  pop_long: pd.DataFrame, 
                                                  gdir: Path):
    """Genera gráficos de incidencia cruda y ASR con IC."""
    from .incidence import incidence_crude_timeseries_with_ci, incidence_asr_timeseries_with_ci
    
    # Cruda
    frames = incidence_crude_timeseries_with_ci(incident, pop_long)
    fig, ax = plt.subplots(figsize=(10.5, 6.4))
    _plot_timeseries_with_ci(ax, frames, "Incidence (crude) by sex and total — 2010–2023", "Rate per 100,000")
    fig.tight_layout()
    fig.savefig(gdir / "incidence_crude_by_sex_total_with_ci.png", dpi=220)
    plt.close(fig)
    
    # WHO ASR
    frames_asr = incidence_asr_timeseries_with_ci(incident, pop_long)
    fig, ax = plt.subplots(figsize=(10.5, 6.4))
    _plot_timeseries_with_ci(ax, frames_asr, "Incidence (WHO-standardized, 80+) by sex and total — 2010–2023", "Rate per 100,000")
    fig.tight_layout()
    fig.savefig(gdir / "incidence_who_asr_by_sex_total_with_ci.png", dpi=220)
    plt.close(fig)


def plot_prevalence_by_year_crude_and_asr_with_ci(base_prev: pd.DataFrame, 
                                                   pop_long: pd.DataFrame, 
                                                   gdir: Path):
    """Genera gráficos de prevalencia cruda y ASR con IC."""
    from .prevalence import prevalence_crude_timeseries_with_ci, prevalence_asr_timeseries_with_ci
    
    # Cruda
    frames = prevalence_crude_timeseries_with_ci(base_prev, pop_long)
    fig, ax = plt.subplots(figsize=(10.5, 6.4))
    _plot_timeseries_with_ci(ax, frames, "Prevalence (crude) by sex and total — 2010–2023", "Rate per 100,000")
    fig.tight_layout()
    fig.savefig(gdir / "prevalence_crude_by_sex_total_with_ci.png", dpi=220)
    plt.close(fig)
    
    # WHO ASR
    frames_asr = prevalence_asr_timeseries_with_ci(base_prev, pop_long)
    fig, ax = plt.subplots(figsize=(10.5, 6.4))
    _plot_timeseries_with_ci(ax, frames_asr, "Prevalence (WHO-standardized, 80+) by sex and total — 2010–2023", "Rate per 100,000")
    fig.tight_layout()
    fig.savefig(gdir / "prevalence_who_asr_by_sex_total_with_ci.png", dpi=220)
    plt.close(fig)


# =============================================================================
# HISTOGRAMA DE PREVALENCIA POR EDAD
# =============================================================================

def prevalence_age_hist_with_totals(base_prev: pd.DataFrame, pop_long: pd.DataFrame, 
                                    year: int, out_path: Path) -> None:
    """
    Histograma de prevalencia por grupo de edad con barras crudas y contribución WHO.
    
    Args:
        base_prev: DataFrame de base de prevalencia
        pop_long: DataFrame de población
        year: Año de referencia
        out_path: Ruta de salida
    """
    alive = base_prev[base_prev["death_is_null"] & (base_prev["diag_year"] <= year)].copy()
    
    if alive.empty:
        fig, ax = plt.subplots(figsize=(13, 6))
        _placeholder(ax, f"Prevalence by age — {year}", "Age group", "Rate per 100,000")
        fig.tight_layout()
        fig.savefig(out_path, dpi=220)
        plt.close(fig)
        return

    cases_age = alive.groupby("edad_grupo").size().rename("cases")
    pop_y = pop_long[pop_long["year"] == year].copy()
    pop_y["edad_grupo"] = pop_y["Edad"].apply(lambda x: to_age_group_5y(x, TOP_OPEN))
    pop_age = pop_y.groupby("edad_grupo")["poblacion"].sum()

    common = sorted(set(cases_age.index) & set(pop_age.index), key=age_group_sortkey)
    
    if not common:
        fig, ax = plt.subplots(figsize=(13, 6))
        _placeholder(ax, f"Prevalence by age — {year}", "Age group", "Rate per 100,000")
        fig.tight_layout()
        fig.savefig(out_path, dpi=220)
        plt.close(fig)
        return

    # Tasa cruda por edad
    crude = pd.Series({
        ag: (cases_age.get(ag, 0) / float(pop_age.get(ag, np.nan)) * 1e5)
        for ag in common
    })

    # Contribución WHO por edad
    who_w = who_world_standard_weights_collapsed()
    weights = pd.Series({ag: who_w[ag] for ag in common if ag in who_w})
    weights = weights / weights.sum() if weights.sum() else weights
    who_contrib = crude[weights.index] * weights

    # Gráfico de barras agrupadas
    fig, ax = plt.subplots(figsize=(13, 6))
    x = np.arange(len(common))
    width = 0.42

    ax.bar(x - width/2, crude.reindex(common).values, width=width, label="Cruda", alpha=0.9)
    ax.bar(x + width/2, who_contrib.reindex(common).values, width=width, label="WHO (contrib.)", alpha=0.9)

    ax.set_xticks(x, common, rotation=45, ha="right")
    ax.set_title(f"Prevalence by age — {year} (Crude vs WHO contribution)", fontsize=14)
    ax.set_xlabel("Age group")
    ax.set_ylabel("Rate per 100,000")
    _poster_axes(ax)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


# =============================================================================
# GRÁFICOS ADICIONALES
# =============================================================================

def plot_new_cases_by_year_by_system(incident: pd.DataFrame, gdir: Path):
    """Genera gráficos de casos nuevos por sistema previsional."""
    years = list(range(INC_START_YEAR, MAX_YEAR + 1))
    inc = incident[incident["year"].isin(years)].copy()
    
    if inc.empty:
        for nm in ["new_cases_by_system_counts.png", "new_cases_by_system_shares.png"]:
            fig, ax = plt.subplots(figsize=(11, 6))
            _placeholder(ax, "New cases by system", "Year", "Count")
            fig.tight_layout()
            fig.savefig(gdir / nm, dpi=220)
            plt.close(fig)
        return
    
    pivot = inc.pivot_table(index="year", columns="prevision", values="sexo", aggfunc="count").fillna(0).astype(int)
    pivot = pivot.reindex(years).fillna(0)
    
    # Conteos
    fig, ax = plt.subplots(figsize=(11, 6))
    pivot.plot(kind="bar", stacked=True, ax=ax)
    ax.set_title("New cases by insurance system — counts", fontsize=14)
    ax.set_xlabel("Year")
    ax.set_ylabel("Count")
    _poster_axes(ax)
    ax.legend(title="System", frameon=False)
    fig.tight_layout()
    fig.savefig(gdir / "new_cases_by_system_counts.png", dpi=220)
    plt.close(fig)
    
    # Shares
    shares = (pivot.div(pivot.sum(axis=1), axis=0) * 100).fillna(0.0)
    fig, ax = plt.subplots(figsize=(11, 6))
    shares.plot(kind="bar", stacked=True, ax=ax)
    ax.set_title("New cases by insurance system — % of yearly total", fontsize=14)
    ax.set_xlabel("Year")
    ax.set_ylabel("% of cases")
    _poster_axes(ax)
    ax.legend(title="System", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    for container in ax.containers:
        ax.bar_label(container, fmt="%.0f%%", label_type="center", color="black")
    fig.tight_layout()
    fig.savefig(gdir / "new_cases_by_system_shares.png", dpi=220)
    plt.close(fig)


def plot_mean_age_at_incidence(incident: pd.DataFrame, gdir: Path):
    """Genera gráfico de edad promedio al diagnóstico por año."""
    years = list(range(INC_START_YEAR, MAX_YEAR + 1))

    ser_total = incident.groupby("year")["edad"].mean().reindex(years)
    ser_f = incident[incident["sexo"] == "F"].groupby("year")["edad"].mean().reindex(years)
    ser_m = incident[incident["sexo"] == "M"].groupby("year")["edad"].mean().reindex(years)

    fig, ax = plt.subplots(figsize=(10, 6))
    has_any = False
    
    for lbl, s in [("Total", ser_total), ("Femenino", ser_f), ("Masculino", ser_m)]:
        if not s.dropna().empty:
            has_any = True
            ax.plot(s.index, s.values, marker="o", linewidth=2.6, label=lbl)
    
    if not has_any:
        _placeholder(ax, "Mean age at incidence by year", "Year", "Mean age (years)")
    else:
        ax.set_title("Mean age at incidence by year (Total / F / M)", fontsize=14)
        ax.set_xlabel("Year")
        ax.set_ylabel("Mean age (years)")
        _poster_axes(ax)
        ax.set_ylim(0, 50)
        ax.set_xticks(years)
        ax.legend(frameon=False, loc="best")
    
    fig.tight_layout()
    fig.savefig(gdir / "mean_age_at_incidence_by_year.png", dpi=220)
    plt.close(fig)


def plot_sex_ratio_and_percent_prevalent(prev_obj: dict, gdir: Path):
    """Genera gráficos de ratio M/F y distribución por sexo."""
    years = list(range(MIN_YEAR, MAX_YEAR + 1))
    counts_sex = prev_obj.get("prev_counts_sex", {})

    # Ratio M/F
    ratio = []
    for y in years:
        s = counts_sex.get(y, pd.Series(dtype=int))
        m = int(s.get("M", 0))
        f = int(s.get("F", 0))
        r = (m / f) if f > 0 else np.nan
        ratio.append(r)
    ratio = pd.Series(ratio, index=years, name="M/F")

    fig, ax = plt.subplots(figsize=(10, 6))
    if ratio.dropna().empty:
        _placeholder(ax, "Sex ratio (M/F) among prevalent", "Year", "Ratio M/F")
    else:
        ax.plot(ratio.index, ratio.values, marker="o", linewidth=2.6)
        ax.set_title("Sex ratio (M/F) among prevalent", fontsize=14)
        ax.set_xlabel("Year")
        ax.set_ylabel("Ratio M/F")
        _poster_axes(ax)
        ax.set_ylim(0, 1)
        ax.set_xticks(years)
    fig.tight_layout()
    fig.savefig(gdir / "sex_ratio_m_over_f_prevalent_by_year.png", dpi=220)
    plt.close(fig)

    # % por sexo
    pctM = []
    pctF = []
    for y in years:
        s = counts_sex.get(y, pd.Series(dtype=int))
        m = int(s.get("M", 0))
        f = int(s.get("F", 0))
        tot = m + f
        pctM.append(m / tot * 100 if tot > 0 else np.nan)
        pctF.append(f / tot * 100 if tot > 0 else np.nan)
    df_pct = pd.DataFrame({"Female": pctF, "Male": pctM}, index=years)

    fig, ax = plt.subplots(figsize=(10, 6))
    if df_pct.dropna().empty:
        _placeholder(ax, "Sex distribution among prevalent — % per year", "Year", "% of prevalent")
    else:
        df_pct.plot(kind="bar", stacked=True, ax=ax)
        ax.set_title("Sex distribution among prevalent — % per year", fontsize=14)
        ax.set_xlabel("Year")
        ax.set_ylabel("% of prevalent")
        _poster_axes(ax)
        ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
        ax.set_xticklabels([str(y) for y in years], rotation=0)
        for container in ax.containers:
            ax.bar_label(container, fmt="%.0f%%", label_type="center", color="black")
    fig.tight_layout()
    fig.savefig(gdir / "sex_percent_prevalent_by_year.png", dpi=220)
    plt.close(fig)


def plot_prevalence_by_system_over_time(base_prev: pd.DataFrame, pop_long: pd.DataFrame, gdir: Path):
    """Genera gráfico de prevalencia por sistema previsional con IC."""
    years = list(range(MIN_YEAR, MAX_YEAR + 1))
    alive = base_prev[base_prev["death_is_null"]].copy()
    systems = ["ISAPRE", "FONASA", "FFAA"]

    rows = []
    for sys in systems:
        alive_sys = alive[alive["prevision"] == sys]
        for y in years:
            d = int((alive_sys["diag_year"] <= y).sum())
            N = float(pop_long[pop_long["year"] == y]["poblacion"].sum())
            r, l, u = poisson_rate_ci(d, N, 1e5)
            rows.append((y, sys, r, l, u))
    df = pd.DataFrame(rows, columns=["year", "system", "rate", "lcl", "ucl"])

    fig, ax = plt.subplots(figsize=(13, 6))
    x = np.arange(len(years))
    width = 0.8 / len(systems)

    for j, sys in enumerate(systems):
        ser = df[df["system"] == sys].set_index("year").reindex(years)
        rates = ser["rate"].to_numpy(dtype=float)
        lcl = ser["lcl"].to_numpy(dtype=float)
        ucl = ser["ucl"].to_numpy(dtype=float)
        xpos = x + (j - (len(systems) - 1) / 2) * width

        ax.bar(xpos, rates, width=width, label=sys, alpha=0.9)
        yerr = np.vstack([np.maximum(0, rates - lcl), np.maximum(0, ucl - rates)])
        ax.errorbar(xpos, rates, yerr=yerr, fmt="none", ecolor="#333333", elinewidth=1.2, capsize=3, capthick=1.2)

    ax.set_xticks(x, years, rotation=45, ha="right")
    ax.set_title("Prevalencia (cruda) por sistema previsional — 2010–2023", fontsize=14)
    ax.set_xlabel("Año")
    ax.set_ylabel("Tasa por 100.000")
    _poster_axes(ax)
    ax.legend(title="Sistema", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(gdir / "prevalence_by_system_over_time_with_ci.png", dpi=220)
    plt.close(fig)

