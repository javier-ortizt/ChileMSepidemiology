"""
Análisis por Sistema Previsional
================================

Funciones para análisis de EM por ISAPRE, FONASA y FFAA:
- Prevalencia cruda y ASR por sistema
- Conteos por ISAPRE específica
- Gráficos detallados por sistema
"""

from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .config import MIN_YEAR, MAX_YEAR, ISAPRE_NAMES, TOP_OPEN
from .utils import normalize, to_age_group_10y
from .statistics import asr_and_ci_from_counts
from .who_standards import who_world_standard_weights_10y
from .plotting import _poster_axes, _placeholder


# =============================================================================
# NORMALIZACIÓN DE NOMBRES DE ISAPRE
# =============================================================================

def _norm_isapre_name(raw: str) -> str | None:
    """
    Normaliza nombre de ISAPRE a nombre estándar.
    
    Args:
        raw: Nombre original de la ISAPRE
        
    Returns:
        Nombre estandarizado o 'ISAPRE (sin detalle)' si no se reconoce
    """
    if raw is None:
        return None
    
    s = normalize(raw)
    
    mapping = {
        "vidatres": "Vida Tres",
        "vida tres": "Vida Tres",
        "cruz blanca": "Cruz Blanca",
        "banmedica": "Banmédica",
        "banmédica": "Banmédica",
        "colmena": "Colmena Golden Cross",
        "colmena golden cross": "Colmena Golden Cross",
        "isapre fundacion": "Isapre Fundación",
        "fundacion": "Isapre Fundación",
        "fundación": "Isapre Fundación",
        "consalud": "Consalud",
        "isalud": "Isalud",
        "nueva masvida": "Nueva Masvida",
        "masvida": "Nueva Masvida",
        "nueva másvida": "Nueva Masvida",
    }
    
    for k, v in mapping.items():
        if k in s:
            return v
    
    return "ISAPRE (sin detalle)"


# =============================================================================
# CONTEOS POR ISAPRE
# =============================================================================

def isapre_counts_by_year(base_prev: pd.DataFrame) -> pd.DataFrame:
    """
    Cuenta prevalentes vivos por ISAPRE específica y año.
    
    Excluye 'ISAPRE (sin detalle)' del resultado.
    
    Args:
        base_prev: DataFrame de base de prevalencia
        
    Returns:
        DataFrame con columnas: year, isapre, count
    """
    if "prevision_raw" not in base_prev.columns:
        return pd.DataFrame(columns=["year", "isapre", "count"])

    alive = base_prev[base_prev["death_is_null"]].copy()
    alive["isapre_name"] = alive["prevision_raw"].apply(_norm_isapre_name)
    alive = alive[alive["isapre_name"].isin(ISAPRE_NAMES)]

    years = range(MIN_YEAR, MAX_YEAR + 1)
    rows = []
    
    for y in years:
        sub = alive[alive["diag_year"] <= y]
        cnt = sub["isapre_name"].value_counts()
        for nm in ISAPRE_NAMES:
            rows.append((y, nm, int(cnt.get(nm, 0))))
    
    return pd.DataFrame(rows, columns=["year", "isapre", "count"])


# =============================================================================
# PREVALENCIA CRUDA POR SISTEMA
# =============================================================================

def crude_prevalence_by_system(base_prev: pd.DataFrame,
                               sys_pop_age: pd.DataFrame,
                               year: int) -> dict:
    """
    Calcula prevalencia cruda por 100,000 por sistema (FONASA/ISAPRE).
    
    Usa como denominador TODA la población del sistema (suma de todos
    los tramos de edad), e incluye en el numerador TODOS los casos
    vivos de ese sistema (incluyendo los sin edad registrada).
    
    Args:
        base_prev: DataFrame de base de prevalencia
        sys_pop_age: DataFrame de población por sistema y edad
        year: Año de referencia
        
    Returns:
        Dict {'FONASA': tasa, 'ISAPRE': tasa} en por 100,000
    """
    # Casos prevalentes vivos hasta ese año
    alive_y = base_prev[base_prev["death_is_null"]].copy()
    alive_y = alive_y[alive_y["diag_year"] <= year].copy()

    # Reclasificación: todo lo que no es FONASA/FFAA -> ISAPRE
    alive_y["prevision2"] = alive_y["prevision"].apply(
        lambda x: "ISAPRE" if x not in ("FONASA", "FFAA") else x
    )

    # Numerador: TODOS los casos vivos por sistema
    num_sys = alive_y["prevision2"].value_counts()

    # Denominador: TODA la población por sistema
    pop_sys = (
        sys_pop_age[sys_pop_age["year"] == year]
        .groupby("prevision")["poblacion"]
        .sum()
    )

    out = {}
    for sys in ["FONASA", "ISAPRE"]:
        cases = int(num_sys.get(sys, 0))
        N = float(pop_sys.get(sys, np.nan))

        if pd.isna(N) or N <= 0:
            out[sys] = np.nan
            print(f"[WARN] Sin denominador válido para {sys} en {year}")
            continue

        rate = cases / N * 1e5
        out[sys] = rate

        print(f"\n=== {sys} — {year} ===")
        print(f"Numerador total (casos vivos)     : {cases}")
        print(f"Denominador total (población)     : {N:,.0f}")
        print(f"Tasa cruda                        : {rate:.6f} por 100.000")

    return out


def debug_crude_prevalence_by_system(base_prev: pd.DataFrame,
                                     sys_pop_age: pd.DataFrame,
                                     year: int) -> pd.DataFrame:
    """
    Versión de debug que muestra numerador y denominador por tramo etario.
    
    Args:
        base_prev: DataFrame de base de prevalencia
        sys_pop_age: DataFrame de población por sistema y edad
        year: Año de referencia
        
    Returns:
        DataFrame con resumen por sistema
    """
    from .utils import age_group_sortkey
    
    alive_y = base_prev[base_prev["death_is_null"]].copy()
    alive_y = alive_y[alive_y["diag_year"] <= year].copy()

    alive_y["prevision2"] = alive_y["prevision"].apply(
        lambda x: "ISAPRE" if x not in ("FONASA", "FFAA") else x
    )
    alive_y["edad_grupo"] = alive_y["edad"].apply(
        lambda x: to_age_group_10y(x, TOP_OPEN)
    )

    rows = []

    for sys in ["FONASA", "ISAPRE"]:
        alive_sys = alive_y[alive_y["prevision2"] == sys]
        cases_age = alive_sys.groupby("edad_grupo").size()

        pop_sys_age = (
            sys_pop_age[
                (sys_pop_age["year"] == year) &
                (sys_pop_age["prevision"] == sys)
            ]
            .set_index("edad_grupo")["poblacion"]
        )

        age_groups = sorted(pop_sys_age.index.unique(), key=age_group_sortkey)

        total_pop = float(pop_sys_age.sum())
        total_cases = int(len(alive_sys))
        rate = total_cases / total_pop * 1e5 if total_pop > 0 else np.nan

        print(f"\n=== {sys} — {year} ===")
        print("Por tramo de edad (numerador / denominador):")
        for ag in age_groups:
            N = float(pop_sys_age.get(ag, 0.0))
            c = int(cases_age.get(ag, 0))
            print(f"  {ag}: {c} / {N:,.0f}")

        print(f"\nTotal casos       : {total_cases}")
        print(f"Total población   : {total_pop:,.0f}")
        print(f"Tasa cruda        : {rate:.6f} por 100.000")

        rows.append({
            "prevision": sys,
            "numerador_casos": total_cases,
            "denominador_poblacion": total_pop,
            "tasa_cruda_por_100k": rate,
        })

    if not rows:
        return pd.DataFrame(
            columns=["numerador_casos", "denominador_poblacion", "tasa_cruda_por_100k"]
        )

    return pd.DataFrame(rows).set_index("prevision")


# =============================================================================
# ASR POR SISTEMA
# =============================================================================

def who_asr_by_system(base_prev: pd.DataFrame,
                      sys_pop_age: pd.DataFrame,
                      year: int) -> dict:
    """
    Calcula prevalencia ajustada por edad (WHO) por sistema con IC95%.
    
    Usa tramos de 10 años para compatibilidad con datos de población.
    
    Args:
        base_prev: DataFrame de base de prevalencia
        sys_pop_age: DataFrame de población por sistema y edad
        year: Año de referencia
        
    Returns:
        Dict {'FONASA': {'ASR': x, 'LCL': x, 'UCL': x}, 'ISAPRE': {...}}
    """
    who_w = who_world_standard_weights_10y()

    alive = base_prev[base_prev["death_is_null"]].copy()
    alive = alive[alive["diag_year"] <= year].copy()
    alive["edad_grupo"] = alive["edad"].apply(lambda x: to_age_group_10y(x, TOP_OPEN))
    alive["prevision2"] = alive["prevision"].apply(
        lambda x: "ISAPRE" if x not in ("FONASA", "FFAA") else x
    )

    results = {}

    for sys in ["FONASA", "ISAPRE"]:
        sub = alive[alive["prevision2"] == sys]
        if sub.empty:
            results[sys] = {"ASR": np.nan, "LCL": np.nan, "UCL": np.nan}
            continue

        cases_age = sub.groupby("edad_grupo").size()

        pop_age = (
            sys_pop_age[
                (sys_pop_age["year"] == year) &
                (sys_pop_age["prevision"] == sys)
            ]
            .set_index("edad_grupo")["poblacion"]
        )

        asr, lcl, ucl = asr_and_ci_from_counts(cases_age, pop_age, who_w)
        results[sys] = {"ASR": asr, "LCL": lcl, "UCL": ucl}

    return results


# =============================================================================
# GRÁFICOS DE ISAPRE
# =============================================================================

def plot_isapre_detail_by_year(base_prev: pd.DataFrame, gdir: Path):
    """
    Genera gráficos de detalle por ISAPRE (conteos y % shares).
    
    Args:
        base_prev: DataFrame de base de prevalencia
        gdir: Directorio de salida de gráficos
    """
    df = isapre_counts_by_year(base_prev)
    
    if df.empty:
        for nm in ["isapre_counts_by_year.png", "isapre_shares_by_year.png"]:
            fig, ax = plt.subplots(figsize=(12, 6))
            _placeholder(ax, "ISAPRE detail by year", "Year", "Count")
            fig.tight_layout()
            fig.savefig(gdir / nm, dpi=220)
            plt.close(fig)
        print("[WARN] No 'prevision_raw' disponible; ISAPRE detail omitido.")
        return

    pivot = (df.pivot_table(index="year", columns="isapre", values="count", aggfunc="sum")
               .reindex(columns=ISAPRE_NAMES).fillna(0).astype(int))

    # Conteos
    fig, ax = plt.subplots(figsize=(12, 6))
    pivot.plot(kind="bar", stacked=True, ax=ax)
    ax.set_title("ISAPRE detail among prevalent (alive) — counts", fontsize=14)
    ax.set_xlabel("Year")
    ax.set_ylabel("Count")
    _poster_axes(ax)
    ax.legend(title="ISAPRE", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.set_xticklabels([str(y) for y in pivot.index], rotation=0)
    fig.tight_layout()
    fig.savefig(gdir / "isapre_counts_by_year.png", dpi=220)
    plt.close(fig)

    # % shares
    shares = (pivot.div(pivot.sum(axis=1), axis=0) * 100).fillna(0.0)
    fig, ax = plt.subplots(figsize=(12, 6))
    shares.plot(kind="bar", stacked=True, ax=ax)
    ax.set_title("ISAPRE detail among prevalent (alive) — % share", fontsize=14)
    ax.set_xlabel("Year")
    ax.set_ylabel("% of prevalent")
    _poster_axes(ax)
    ax.legend(title="ISAPRE", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    ax.set_xticklabels([str(y) for y in shares.index], rotation=0)
    for container in ax.containers:
        ax.bar_label(container, fmt="%.0f%%", label_type="center", color="black")
    fig.tight_layout()
    fig.savefig(gdir / "isapre_shares_by_year.png", dpi=220)
    plt.close(fig)

