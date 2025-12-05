"""
Análisis de Correlación
=======================

Correlación entre prevalencia y latitud geográfica.
"""

from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

from .config import MIN_YEAR, MAX_YEAR
from .statistics import pearson_ci, ols_confidence_band, asr_from_rates_by_group
from .who_standards import who_world_standard_weights_collapsed
from .geography import _region_lookup_from_pop, get_latitudes_df
from .prevalence import prevalence_by_region_age_for_asr
from .plotting import _placeholder, _poster_axes


# =============================================================================
# CORRELACIÓN PEARSON: PREVALENCIA VS LATITUD
# =============================================================================

def pearson_prevalence_vs_latitude(prev_obj: dict,
                                   pop_long: pd.DataFrame,
                                   base_prev: pd.DataFrame,
                                   year: int,
                                   outdir: Path,
                                   gdir: Path) -> None:
    """
    Calcula y grafica correlación de Pearson entre prevalencia y latitud.
    
    Genera dos gráficos:
    - Prevalencia cruda vs latitud
    - Prevalencia ASR (WHO) vs latitud
    
    Incluye banda de confianza 95% y estadísticas en el gráfico.
    
    Args:
        prev_obj: Dict de conteos de prevalencia
        pop_long: DataFrame de población INE
        base_prev: DataFrame de base de prevalencia
        year: Año de referencia
        outdir: Directorio para CSV de salida
        gdir: Directorio para gráficos
    """
    # Verificar datos disponibles
    if not prev_obj.get("prev_counts_region"):
        for name, ylab in [
            (f"prevalence_pearson_crude_vs_latitude_{year}.png", "Prevalence per 100,000"),
            (f"prevalence_pearson_who_asr_vs_latitude_{year}.png", "WHO-standardized prevalence (per 100,000)")
        ]:
            fig, ax = plt.subplots(figsize=(8.2, 6))
            _placeholder(ax, f"No crude/ASR data for {year}", "Latitude (°S)", ylab)
            fig.tight_layout()
            fig.savefig(gdir / name, dpi=220)
            plt.close(fig)
        print(f"\n[INFO] No regional counts found for {year}. Placeholder plots saved.")
        return

    counts_series = prev_obj["prev_counts_region"].get(year, pd.Series(dtype=int))
    if counts_series is None or counts_series.empty:
        for name, ylab in [
            (f"prevalence_pearson_crude_vs_latitude_{year}.png", "Prevalence per 100,000"),
            (f"prevalence_pearson_who_asr_vs_latitude_{year}.png", "WHO-standardized prevalence (per 100,000)")
        ]:
            fig, ax = plt.subplots(figsize=(8.2, 6))
            _placeholder(ax, f"No regional data for {year}", "Latitude (°S)", ylab)
            fig.tight_layout()
            fig.savefig(gdir / name, dpi=220)
            plt.close(fig)
        print(f"\n[INFO] No regional prevalence for {year}. Placeholder plots saved.")
        return

    counts_series.index.name = "region_code"
    counts_series = counts_series.rename("cases")

    # Calcular prevalencia cruda por región
    pop_y = pop_long[pop_long["year"] == year].groupby("region_code")["poblacion"].sum()
    common_idx = counts_series.index.intersection(pop_y.index)
    crude_rate = (counts_series.loc[common_idx] / pop_y.loc[common_idx] * 1e5)
    crude_rate.name = "crude_prev_per_100k"

    # Calcular WHO ASR por región
    asr_reg = pd.Series(dtype=float, name="who_asr_per_100k")
    reg_age_df = prevalence_by_region_age_for_asr(base_prev, pop_long)
    
    if not reg_age_df.empty:
        asr_all = asr_from_rates_by_group(reg_age_df, who_world_standard_weights_collapsed(), "region_code")
        if not asr_all.empty:
            try:
                asr_year = asr_all.xs(year, level=0)
                asr_year.index.name = "region_code"
                asr_reg = asr_year.rename("who_asr_per_100k")
            except KeyError:
                pass

    # Unir con nombres de región y latitudes
    look = _region_lookup_from_pop(pop_long)
    df_rates = look.merge(crude_rate.rename("crude_prev_per_100k"),
                          left_on="region_code", right_index=True, how="left")
    if not asr_reg.empty:
        df_rates = df_rates.merge(asr_reg.rename("who_asr_per_100k"),
                                  left_on="region_code", right_index=True, how="left")
    
    lat_tbl = get_latitudes_df()[["REGION_NORM", "latitud"]]
    merged = df_rates.merge(lat_tbl, on="REGION_NORM", how="left")

    # Exportar CSV
    out_csv = outdir / f"prevalence_latitude_merge_{year}.csv"
    merged.to_csv(out_csv, index=False, encoding="utf-8-sig")

    # Función de graficado
    def _scatter_band(ax, x, y, title, ylab):
        r, p = pearsonr(x, y)
        ci_low, ci_high = pearson_ci(r, len(x))
        
        ax.scatter(x, y, s=48, alpha=0.9)
        
        xfit = np.linspace(x.min(), x.max(), 200)
        yfit, lcl, ucl = ols_confidence_band(x, y, xfit)
        
        ax.plot(xfit, yfit, linewidth=2.2)
        ax.fill_between(xfit, lcl, ucl, alpha=0.15, linewidth=0)
        
        ax.set_title(title, fontsize=12)
        ax.set_xlabel("Latitude (°S)")
        ax.set_ylabel(ylab)
        _poster_axes(ax)
        
        # Anotación en el gráfico
        txt = f"r = {r:.3f}\n95% CI = [{ci_low:.3f}, {ci_high:.3f}]\np = {p:.3g}\nn = {len(x)}"
        ax.text(0.98, 0.02, txt, transform=ax.transAxes, ha="right", va="bottom",
                fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#999999", alpha=0.9))

        print(f"{title}: r={r:.3f} (95% CI {ci_low:.3f}–{ci_high:.3f}), p={p:.4g}, n={len(x)}")

    # Gráfico: Prevalencia CRUDA vs latitud
    fig, ax = plt.subplots(figsize=(8.6, 6.2))
    m1 = merged.dropna(subset=["crude_prev_per_100k", "latitud"])
    
    if len(m1) >= 3:
        x = m1["latitud"].astype(float).to_numpy()
        y = m1["crude_prev_per_100k"].astype(float).to_numpy()
        _scatter_band(ax, x, y, f"Crude prevalence vs latitude ({year})", "Prevalence per 100,000")
    else:
        _placeholder(ax, f"Crude prevalence vs latitude ({year})", "Latitude (°S)", "Prevalence per 100,000")
        print(f"\n[INFO] Insufficient data for crude Pearson {year} (n={len(m1)})")
    
    fig.tight_layout()
    fig.savefig(gdir / f"prevalence_pearson_crude_vs_latitude_{year}.png", dpi=220)
    plt.close(fig)

    # Gráfico: Prevalencia WHO ASR vs latitud
    fig, ax = plt.subplots(figsize=(8.6, 6.2))
    
    if "who_asr_per_100k" in merged.columns:
        m2 = merged.dropna(subset=["who_asr_per_100k", "latitud"])
    else:
        m2 = pd.DataFrame(columns=["who_asr_per_100k", "latitud"])
    
    if len(m2) >= 3:
        x = m2["latitud"].astype(float).to_numpy()
        y = m2["who_asr_per_100k"].astype(float).to_numpy()
        _scatter_band(ax, x, y, f"WHO-standardized prevalence vs latitude ({year})", 
                     "WHO-standardized prevalence (per 100,000)")
    else:
        _placeholder(ax, f"WHO-standardized prevalence vs latitude ({year})",
                     "Latitude (°S)", "WHO-standardized prevalence (per 100,000)")
        print(f"[INFO] Insufficient data for WHO ASR Pearson {year} (n={len(m2)})")
    
    fig.tight_layout()
    fig.savefig(gdir / f"prevalence_pearson_who_asr_vs_latitude_{year}.png", dpi=220)
    plt.close(fig)

    print(f"(Merge saved to: {out_csv})")

