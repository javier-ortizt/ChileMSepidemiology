"""
Carga de Datos (I/O)
====================

Funciones para leer los archivos de entrada:
- Base de datos de EM
- Proyecciones poblacionales INE
- Población por sistema previsional
"""

from __future__ import annotations
from pathlib import Path
import re
import warnings

import pandas as pd
import numpy as np

from .config import (
    EXCEL_FILE, INE_POP_FILE, PREV_REGION_FILE, BEN_SYS_FILE,
    MAX_YEAR, TOP_OPEN
)
from .utils import normalize, age_group_sortkey, to_age_group_10y
from .geography import normalizar_nombre_region, _fix_region_norm_aliases, _region_lookup_from_pop


# =============================================================================
# CARGA DE BASE DE DATOS PRINCIPAL
# =============================================================================

def read_bdem(filepath: str = EXCEL_FILE) -> pd.DataFrame:
    """
    Lee el archivo Excel de la base de datos de EM.
    
    Args:
        filepath: Ruta al archivo Excel
        
    Returns:
        DataFrame con los datos de casos
        
    Raises:
        FileNotFoundError: Si el archivo no existe
    """
    if not Path(filepath).exists():
        raise FileNotFoundError(f"Archivo no encontrado: {filepath}")
    return pd.read_excel(filepath, engine="openpyxl")


# =============================================================================
# CARGA DE PROYECCIONES POBLACIONALES INE
# =============================================================================

def load_ine_population(path: str = INE_POP_FILE) -> pd.DataFrame:
    """
    Carga las proyecciones poblacionales del INE (2002-2035) por comuna.
    
    Retorna formato largo con columnas:
    - year: Año
    - sexo: M/F
    - Edad: Edad simple
    - region_code: Código de región (int)
    - region_name: Nombre de región
    - poblacion: Población
    
    Args:
        path: Ruta al archivo de proyecciones
        
    Returns:
        DataFrame en formato largo
        
    Raises:
        FileNotFoundError: Si el archivo no existe
        ValueError: Si faltan columnas requeridas
    """
    if not Path(path).exists():
        raise FileNotFoundError(f"Archivo de población INE no encontrado: {path}")
    
    df = pd.read_excel(path, engine="openpyxl")

    # Identificar columna de sexo -> M/F
    sexo_col = next((c for c in df.columns if "sexo" in normalize(c)), None)
    df["__sexo"] = df[sexo_col].map({1: "M", 2: "F"}).fillna(df[sexo_col])

    # Identificar código y nombre de región
    region_code_col = next((c for c in df.columns if normalize(c) == "region"), None)
    region_name_col = next((c for c in df.columns if normalize(c) in ("nombre region", "nombreregion")), None)
    
    if region_code_col is None:
        raise ValueError("Columna 'Region' (código) no encontrada en archivo INE.")
    
    df["__region_code"] = pd.to_numeric(df[region_code_col], errors="coerce").astype("Int64")
    df["__region_name"] = df[region_name_col].astype(str) if region_name_col else pd.Series(pd.NA, index=df.index)

    # Identificar columnas de población por año
    pop_cols = [c for c in df.columns if normalize(c).startswith("poblacion ")]
    
    # Convertir a formato largo
    long = df.melt(
        id_vars=["__sexo", "Edad", "__region_code", "__region_name"],
        value_vars=pop_cols,
        var_name="__year_col",
        value_name="poblacion"
    )
    
    # Extraer año de nombre de columna
    years = long["__year_col"].astype(str).str.extract(r"(\d{4})")[0]
    long = long.assign(year=years).dropna(subset=["year"]).copy()
    long["year"] = long["year"].astype(int)

    # Agregar por (año, región, sexo, edad)
    long = long.groupby(
        ["year", "__region_code", "__region_name", "__sexo", "Edad"],
        as_index=False
    )["poblacion"].sum()
    
    long.rename(columns={
        "__sexo": "sexo",
        "__region_code": "region_code",
        "__region_name": "region_name"
    }, inplace=True)

    return long


# =============================================================================
# CARGA DE POBLACIÓN POR PREVISIÓN Y REGIÓN
# =============================================================================

def load_population_by_prevision_region(path: str = PREV_REGION_FILE,
                                        pop_long: pd.DataFrame | None = None,
                                        year: int = MAX_YEAR) -> pd.DataFrame:
    """
    Lee población por previsión y región desde 'Previsionporregion.xlsx'.
    
    Retorna DataFrame largo con:
    - year: Año
    - region_code: Código de región (int)
    - prevision: 'FONASA' o 'ISAPRE'
    - poblacion: Población
    
    Args:
        path: Ruta al archivo
        pop_long: DataFrame de población INE (para mapear región)
        year: Año de referencia
        
    Returns:
        DataFrame con población por previsión y región
    """
    if not Path(path).exists():
        raise FileNotFoundError(f"Población por previsión no encontrada: {path}")

    df = pd.read_excel(path, sheet_name=0, engine="openpyxl")
    df = df.rename(columns={c: str(c).strip() for c in df.columns})
    
    if "Region" not in df.columns:
        raise ValueError("Se esperaba columna 'Region' en Previsionporregion.xlsx")

    # Identificar columnas de sistemas
    value_cols = [c for c in df.columns if c.upper() in ("FONASA", "ISAPRE", "FFAA")]
    if not value_cols:
        raise ValueError("Previsionporregion.xlsx no contiene columnas de población (FONASA/ISAPRE/FFAA).")

    # Convertir a formato largo
    long = df.melt(
        id_vars=["Region"],
        value_vars=value_cols,
        var_name="prevision",
        value_name="poblacion"
    ).dropna(subset=["poblacion"])
    
    long["prevision"] = long["prevision"].str.upper().str.strip()
    long["poblacion"] = pd.to_numeric(long["poblacion"], errors="coerce").fillna(0).astype(int)

    # Mapear Region -> region_code usando INE
    if pop_long is None or pop_long.empty:
        raise ValueError("Se requiere 'pop_long' (INE) para mapear Region -> region_code.")

    look = _region_lookup_from_pop(pop_long)[["region_code", "REGION_NORM"]]
    tmp = long.copy()
    tmp["REGION_NORM"] = tmp["Region"].map(normalizar_nombre_region).pipe(_fix_region_norm_aliases)
    merged = tmp.merge(look, on="REGION_NORM", how="left")

    # Validar mapeo
    miss = merged[merged["region_code"].isna()]
    if not miss.empty:
        warnings.warn(f"[WARN] No se pudo mapear algunas regiones:\n{miss['Region'].unique()}")

    merged = merged.dropna(subset=["region_code"]).copy()
    merged["region_code"] = merged["region_code"].astype(int)
    merged["year"] = int(year)

    return merged[["year", "region_code", "prevision", "poblacion"]].copy()


# =============================================================================
# CARGA DE POBLACIÓN POR SISTEMA Y TRAMO ETARIO
# =============================================================================

def load_population_by_system_age(path: str = BEN_SYS_FILE,
                                  year: int = MAX_YEAR) -> pd.DataFrame:
    """
    Lee población de beneficiarios por sistema y tramo etario.
    
    Args:
        path: Ruta al archivo
        year: Año de referencia
        
    Returns:
        DataFrame con columnas: year, edad_grupo, prevision, poblacion
    """
    df = pd.read_excel(path, engine="openpyxl")
    df = df.rename(columns=lambda c: str(c).strip())

    age_col = "TRAMOETARIO"

    def _norm_tramo(s: str) -> str:
        """Normaliza tramo etario a formato estándar."""
        s = str(s).lower()
        s = s.replace("años", "").replace("anos", "")
        s = s.replace("y más", "+").replace("y mas", "+").replace("ymas", "+")
        s = s.replace(" a ", "-").replace("a ", "-").replace(" a", "-")
        s = re.sub(r"\s+", "", s)
    
        # Si no hay dígitos, es desconocido
        if not re.search(r"\d", s):
            return "Unknown"
    
        if "+" in s:
            return f"{TOP_OPEN}+"
    
        m = re.match(r"(\d+)-(\d+)", s)
        if m:
            lo, hi = int(m.group(1)), int(m.group(2))
            if hi >= TOP_OPEN:
                return f"{TOP_OPEN}+"
            return f"{lo}-{hi}"
    
        # Número suelto -> grupo de 10 años
        try:
            a = int(re.findall(r"\d+", s)[0])
            return to_age_group_10y(a, TOP_OPEN)
        except Exception:
            return "Unknown"

    df["edad_grupo"] = df[age_col].astype(str).apply(_norm_tramo)

    # Identificar columnas de sistemas
    sys_cols = [c for c in df.columns if c.upper() in ("FONASA", "ISAPRE")]
    if not sys_cols:
        raise ValueError("No se encontraron columnas FONASA/ISAPRE en beneficiariosTramoEtario.")

    # Convertir a formato largo
    long = df.melt(
        id_vars=["edad_grupo"],
        value_vars=sys_cols,
        var_name="prevision",
        value_name="poblacion",
    )

    long["prevision"] = long["prevision"].str.upper()
    long["poblacion"] = pd.to_numeric(long["poblacion"], errors="coerce").fillna(0)
    long["year"] = int(year)

    return long[["year", "edad_grupo", "prevision", "poblacion"]]

