"""
Mapas Coropléticos
==================

Funciones para generar mapas de prevalencia por región.
"""

from __future__ import annotations
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import FuncFormatter
import geopandas as gpd

from .config import SHAPEFILE_PATH, MAX_YEAR
from .utils import fmt_lon, fmt_lat
from .geography import (
    normalizar_nombre_region, _fix_region_norm_aliases,
    _shape_region_code_col, _shape_region_name_col,
    _region_name_table_from_pop
)


# =============================================================================
# CÁLCULO DE PREVALENCIA REGIONAL
# =============================================================================

def _compute_regional_crude_prev(base_prev: pd.DataFrame,
                                 pop_long: pd.DataFrame,
                                 year: int,
                                 filter_prevision: str | None = None,
                                 sys_pop_df: pd.DataFrame | None = None) -> pd.Series:
    """
    Calcula prevalencia cruda por 100,000 por región para un año dado.
    
    Args:
        base_prev: DataFrame de base de prevalencia
        pop_long: DataFrame de población INE
        year: Año de referencia
        filter_prevision: Filtrar por sistema previsional (FONASA/ISAPRE)
        sys_pop_df: DataFrame de población por sistema y región (denominador)
        
    Returns:
        Serie con región como índice y tasa como valor
    """
    if base_prev is None or base_prev.empty:
        return pd.Series(dtype=float)

    alive = base_prev[base_prev["death_is_null"]].copy()
    alive = alive[(alive["diag_year"].notna()) & (alive["diag_year"] <= year)]
    alive = alive.dropna(subset=["region_code"])
    alive["region_code"] = alive["region_code"].astype(int)

    # Filtrar por previsión si corresponde
    if filter_prevision:
        alive = alive[alive["prevision"].str.upper() == str(filter_prevision).upper()]

    cases_reg = alive.groupby("region_code").size().rename("cases")

    # Denominador
    if filter_prevision and sys_pop_df is not None:
        sys = str(filter_prevision).upper()
        pop_y = (sys_pop_df[(sys_pop_df["year"] == year) & (sys_pop_df["prevision"] == sys)]
                 .groupby("region_code")["poblacion"].sum())
        
        if pop_y.empty and sys != "FFAA":
            warnings.warn(f"[WARN] No hay población para {sys} en {year}; se omite mapa.")
            return pd.Series(dtype=float)
        
        if pop_y.empty and sys == "FFAA":
            warnings.warn("[WARN] 'FFAA' no disponible; se usa población general como aproximación.")
            pop_y = (pop_long[pop_long["year"] == year]
                     .dropna(subset=["region_code"])
                     .assign(region_code=lambda d: d["region_code"].astype(int))
                     .groupby("region_code")["poblacion"].sum())
    else:
        pop_y = (pop_long[pop_long["year"] == year]
                 .dropna(subset=["region_code"])
                 .assign(region_code=lambda d: d["region_code"].astype(int))
                 .groupby("region_code")["poblacion"].sum())

    common = cases_reg.index.intersection(pop_y.index)
    if len(common) == 0:
        return pd.Series(dtype=float)

    rates = (cases_reg.loc[common] / pop_y.loc[common] * 1e5).astype(float)
    rates.name = "prev_per_100k"
    return rates


# =============================================================================
# MAPA COROPLÉTICO
# =============================================================================

def plot_prevalence_map_by_region(base_prev: pd.DataFrame,
                                  pop_long: pd.DataFrame,
                                  shapefile_path: str,
                                  year: int,
                                  out_path: Path,
                                  title: str | None = None,
                                  filter_prevision: str | None = None,
                                  cmap_name: str = "Reds",
                                  sys_pop_df: pd.DataFrame | None = None) -> None:
    """
    Genera mapa coroplético de prevalencia cruda por región.
    
    Args:
        base_prev: DataFrame de base de prevalencia
        pop_long: DataFrame de población INE
        shapefile_path: Ruta al shapefile de regiones
        year: Año de referencia
        out_path: Ruta de salida del gráfico
        title: Título del mapa (opcional)
        filter_prevision: Filtrar por sistema (FONASA/ISAPRE)
        cmap_name: Nombre del colormap (Reds, Blues, Greens)
        sys_pop_df: DataFrame de población por sistema y región
    """
    rates = _compute_regional_crude_prev(base_prev, pop_long, year,
                                         filter_prevision=filter_prevision,
                                         sys_pop_df=sys_pop_df)
    
    if rates.empty:
        print(f"[WARN] No regional data for map {title or ''}. Saving placeholder.")
        fig, ax = plt.subplots(1, 1, figsize=(9, 11))
        ax.text(0.5, 0.5, "No regional data", ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
        fig.tight_layout()
        fig.savefig(out_path, dpi=220)
        plt.close(fig)
        return

    # Cargar shapefile
    gdf = gpd.read_file(shapefile_path).to_crs(epsg=4326)

    # Unir datos por código o nombre de región
    code_col = _shape_region_code_col(gdf)
    if code_col:
        gdf["_region_code_"] = pd.to_numeric(gdf[code_col], errors="coerce").astype("Int64")
        merged = gdf.merge(rates.rename("rate"), left_on="_region_code_", right_index=True, how="left")
    else:
        name_col = _shape_region_name_col(gdf)
        gdf["REGION_NORM"] = gdf[name_col].map(normalizar_nombre_region).pipe(_fix_region_norm_aliases)
        nm = _region_name_table_from_pop(pop_long)
        nm["REGION_NORM"] = nm["region_name"].map(normalizar_nombre_region).pipe(_fix_region_norm_aliases)
        df_join = nm.merge(rates.rename("rate"), left_on="region_code", right_index=True, how="right")
        merged = gdf.merge(df_join[["REGION_NORM", "rate"]], on="REGION_NORM", how="left")

    vals = merged["rate"].astype(float)
    if vals.dropna().empty:
        print("[WARN] All regional rates are NaN after merge; map skipped.")
        return

    # Configurar escala de colores
    vmin, vmax = 0.0, float(np.nanmax(vals))
    boundaries = np.linspace(vmin, vmax, 6)
    norm = BoundaryNorm(boundaries, ncolors=256)

    # Crear figura
    fig, ax = plt.subplots(1, 1, figsize=(5, 11))
    merged.plot(column="rate", cmap=cmap_name, linewidth=0.6,
                edgecolor="black", ax=ax, legend=False, norm=norm)
    merged.boundary.plot(ax=ax, linewidth=0.8, edgecolor="black")
    
    # Límites geográficos de Chile
    ax.set_xlim(-80.2, -65.0)
    ax.set_ylim(-58.0, -15.5)

    # Formateo de ejes
    ax.xaxis.set_major_formatter(FuncFormatter(fmt_lon))
    ax.yaxis.set_major_formatter(FuncFormatter(fmt_lat))
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(True, linestyle="--", alpha=0.4)

    # Título
    ttl = title or f"Crude prevalence by region (per 100,000) — {year}"
    if filter_prevision:
        ttl += f" — {filter_prevision.title()}"
    ax.set_title(ttl, fontsize=14)

    # Barra de colores
    smap = plt.cm.ScalarMappable(cmap=cmap_name, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    smap._A = []
    cbar = fig.colorbar(smap, ax=ax, orientation='vertical', fraction=0.025, pad=0.02, shrink=0.85)
    cbar.set_label('Prevalence (per 100,000)')
    cbar.set_ticks(boundaries)
    cbar.ax.set_yticklabels([f"{x:.0f}" for x in boundaries])
    
    fig.subplots_adjust(top=0.90, bottom=0.05, left=0.05, right=0.88)
    fig.tight_layout(rect=[0, 0, 0.88, 0.93])

    fig.savefig(out_path, dpi=220)
    plt.close(fig)

    print(f"[OK] Map saved: {out_path}")

