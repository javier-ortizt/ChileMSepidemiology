"""
Helpers Geográficos
===================

Funciones para manejo de regiones, normalización de nombres
y coordenadas de latitud para análisis geográfico.
"""

from __future__ import annotations
import re

import pandas as pd
import unidecode


# =============================================================================
# NORMALIZACIÓN DE NOMBRES DE REGIÓN
# =============================================================================

def normalizar_nombre_region(nombre: str) -> str:
    """
    Normaliza nombre de región para comparaciones robustas.
    
    Args:
        nombre: Nombre de región original
        
    Returns:
        Nombre normalizado (sin acentos, minúsculas, sin palabras auxiliares)
    """
    s = str(nombre).lower().strip()
    s = unidecode.unidecode(s)
    for w in ["region", "región", "del", "de", ".", ","]:
        s = s.replace(w, " ")
    s = re.sub(r"\s+", " ", s).strip()
    return s


def _fix_region_norm_aliases(s: pd.Series) -> pd.Series:
    """
    Corrige alias conocidos de nombres de región.
    
    Args:
        s: Serie con nombres normalizados
        
    Returns:
        Serie con alias corregidos
    """
    repl = {
        "biobio": "bio-bio",
        "libertador b o'higgins": "libertador bernardo o'higgins",
        "magallanes y la antartica chilena": "magallanes y antartica chilena",
        "aisen gral c ibanez campo": "aysen gralibanez campo",
    }
    return s.replace(repl)


# =============================================================================
# IDENTIFICACIÓN DE COLUMNAS EN SHAPEFILES
# =============================================================================

def _shape_region_code_col(gdf) -> str | None:
    """
    Encuentra la columna de código de región en un GeoDataFrame.
    
    Args:
        gdf: GeoDataFrame del shapefile
        
    Returns:
        Nombre de la columna o None
    """
    from .utils import normalize
    
    cands = ["cod_reg", "cod region", "region", "codigo region", "codregion", "cod_deis_region", "cod_deis"]
    normcols = {normalize(c): c for c in gdf.columns}
    
    for want in cands:
        for ncol, orig in normcols.items():
            if want in ncol:
                try:
                    pd.to_numeric(gdf[orig], errors="coerce")
                    return orig
                except Exception:
                    continue
    return None


def _shape_region_name_col(gdf) -> str | None:
    """
    Encuentra la columna de nombre de región en un GeoDataFrame.
    
    Args:
        gdf: GeoDataFrame del shapefile
        
    Returns:
        Nombre de la columna o None
    """
    from .utils import normalize
    
    cands = ["region", "nom_reg", "nombre region", "nombre_region", "nom region"]
    normcols = {normalize(c): c for c in gdf.columns}
    
    for want in cands:
        for ncol, orig in normcols.items():
            if want in ncol:
                return orig
    
    for c in gdf.columns:
        if str(c).lower() == "region":
            return c
    return None


# =============================================================================
# TABLAS DE LOOKUP REGIÓN <-> CÓDIGO
# =============================================================================

def _region_name_table_from_pop(pop_long: pd.DataFrame) -> pd.DataFrame:
    """
    Construye tabla de mapeo region_code -> region_name desde población INE.
    
    Args:
        pop_long: DataFrame de población INE
        
    Returns:
        DataFrame con region_code, region_name, region_name_norm
    """
    df = pop_long.dropna(subset=["region_code"]).copy()
    df["region_code"] = df["region_code"].astype(int)
    
    nm = (df.groupby(["region_code", "region_name"]).size()
            .reset_index(name="n")
            .sort_values(["region_code", "n"], ascending=[True, False]))
    nm = nm.drop_duplicates("region_code")[["region_code", "region_name"]]
    nm["region_name_norm"] = nm["region_name"].map(normalizar_nombre_region)
    
    return nm


def _region_lookup_from_pop(pop_long: pd.DataFrame) -> pd.DataFrame:
    """
    Construye tabla de lookup con nombres normalizados y corregidos.
    
    Args:
        pop_long: DataFrame de población INE
        
    Returns:
        DataFrame con region_code, region_name, REGION_NORM
    """
    df = pop_long.dropna(subset=["region_code"]).copy()
    df["region_code"] = df["region_code"].astype(int)
    
    nm = (df.groupby(["region_code", "region_name"]).size()
            .reset_index(name="n")
            .sort_values(["region_code", "n"], ascending=[True, False]))
    nm = nm.drop_duplicates("region_code")[["region_code", "region_name"]]
    nm["REGION_NORM"] = nm["region_name"].map(normalizar_nombre_region).pipe(_fix_region_norm_aliases)
    
    return nm


# =============================================================================
# TABLA DE LATITUDES POR REGIÓN
# =============================================================================

def get_latitudes_df() -> pd.DataFrame:
    """
    Retorna latitudes aproximadas del centroide de cada región.
    
    Returns:
        DataFrame con REGION_NORM y latitud (grados Sur, valores positivos)
    """
    data = [
        ("Arica y Parinacota",                          18.5),
        ("Tarapacá",                                    20.2),
        ("Antofagasta",                                 23.6),
        ("Atacama",                                     27.4),
        ("Coquimbo",                                    30.0),
        ("Valparaíso",                                  32.9),
        ("Metropolitana de Santiago",                   33.5),
        ("Libertador Bernardo O'Higgins",               34.6),
        ("Maule",                                       35.5),
        ("Ñuble",                                       36.7),
        ("Biobío",                                      37.4),
        ("La Araucanía",                                38.5),
        ("Los Ríos",                                    40.0),
        ("Los Lagos",                                   41.5),
        ("Aysén del Gral. Carlos Ibáñez del Campo",     45.5),
        ("Magallanes y Antártica Chilena",              53.2),
    ]
    
    df = pd.DataFrame(data, columns=["region_name", "latitud"])
    df["REGION_NORM"] = df["region_name"].map(normalizar_nombre_region).pipe(_fix_region_norm_aliases)
    
    return df[["REGION_NORM", "latitud"]]

