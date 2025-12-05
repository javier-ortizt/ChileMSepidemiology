"""
Cálculo de Incidencia
=====================

Funciones para el cálculo de tasas de incidencia de EM:
- Extracción de casos incidentes (primer año por beneficiario)
- Tasas crudas y específicas por edad
- Series temporales con intervalos de confianza
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .config import INC_START_YEAR, MAX_YEAR, TOP_OPEN
from .utils import (
    find_column, classify_prevision, normalize_sex, 
    to_age_group_5y, sort_age_index
)
from .statistics import poisson_rate_ci, asr_and_ci_from_counts
from .who_standards import who_world_standard_weights_collapsed


# =============================================================================
# EXTRACCIÓN DE AÑO DEL CASO
# =============================================================================

def extract_year_column(df: pd.DataFrame) -> pd.Series | None:
    """
    Extrae el año de creación del caso GES desde la columna correspondiente.
    
    Busca primero columna de año directo, luego columna de fecha.
    
    Args:
        df: DataFrame con datos de casos
        
    Returns:
        Serie con años o None si no se encuentra
    """
    # Buscar columna de año
    col_year = find_column(df, ["año creacion caso ges", "ano creacion caso ges", "anio creacion caso ges"])
    if col_year and col_year in df.columns:
        return pd.to_numeric(df[col_year], errors="coerce").astype("Int64")
    
    # Buscar columna de fecha
    col_date = find_column(df, ["fecha creacion caso ges", "fecha creacion"])
    if col_date and col_date in df.columns:
        dates = pd.to_datetime(df[col_date], errors="coerce", dayfirst=True)
        return dates.dt.year.astype("Int64")
    
    return None


# =============================================================================
# CASOS INCIDENTES
# =============================================================================

def compute_incidence_from_bdem(df_cases: pd.DataFrame) -> pd.DataFrame:
    """
    Identifica casos incidentes (primer año por beneficiario) con ajustes.
    
    Ajustes aplicados:
    - FFAA con año < 2010 se asignan a 2010 (baseline)
    - 2010 se excluye del reporte de incidencia (es baseline)
    
    Args:
        df_cases: DataFrame con todos los registros de casos
        
    Returns:
        DataFrame con casos incidentes únicos:
        - year: Año de incidencia
        - sexo: M/F
        - edad: Edad al diagnóstico
        - edad_grupo: Grupo quinquenal
        - region_code: Código de región
        - prevision: FONASA/ISAPRE/FFAA
        
    Raises:
        ValueError: Si no se puede determinar el año del caso
    """
    year_series = extract_year_column(df_cases)
    if year_series is None:
        raise ValueError("No se pudo determinar el año del caso.")
    
    dfy = df_cases.copy()
    dfy["__year"] = year_series

    # Clasificar previsión
    prev_col = find_column(dfy, ["prevision"])
    if prev_col and prev_col in dfy.columns:
        dfy["__prev"] = dfy[prev_col].astype(str).apply(classify_prevision)
    else:
        dfy["__prev"] = np.nan

    # Ajuste FFAA: años 2001-2009 -> 2010
    dfy["__year_adj"] = dfy["__year"]
    mask_ffaa_pre2010 = (dfy["__prev"] == "FFAA") & dfy["__year"].notna() & (dfy["__year"] < 2010)
    dfy.loc[mask_ffaa_pre2010, "__year_adj"] = 2010

    # Identificar columnas relevantes
    ben_id = find_column(dfy, ["cod beneficiario", "id beneficiario", "id", "codigo beneficiario", "rut", "n°", "n "])
    col_sexo = find_column(dfy, ["sexo beneficiario", "sexo"])
    col_edad = find_column(dfy, ["edad (anos) creacion caso ges", "edad (años) creacion caso ges", "edad"])
    region_code_col = find_column(dfy, ["cod. deis region beneficiario", "cod deis region beneficiario", "cod deis region"])

    cols = [ben_id, "__year_adj", col_sexo, col_edad, "__prev"]
    if region_code_col:
        cols.append(region_code_col)

    # Filtrar y ordenar
    tmp = dfy[cols].dropna(subset=[ben_id, "__year_adj"]).copy()
    tmp.sort_values([ben_id, "__year_adj"], inplace=True)
    
    # Mantener solo primera aparición de cada beneficiario
    incident = tmp.drop_duplicates(subset=[ben_id], keep="first").copy()

    # Renombrar columnas
    rename_map = {
        "__year_adj": "year",
        col_sexo: "sexo",
        col_edad: "edad",
        "__prev": "prevision"
    }
    if region_code_col:
        rename_map[region_code_col] = "region_code"
    incident.rename(columns=rename_map, inplace=True)

    # Normalizar valores
    incident["sexo"] = incident["sexo"].apply(normalize_sex)
    incident["edad_grupo"] = incident["edad"].apply(to_age_group_5y)
    incident["year"] = pd.to_numeric(incident["year"], errors="coerce")
    
    if "region_code" in incident.columns:
        incident["region_code"] = pd.to_numeric(incident["region_code"], errors="coerce").astype("Int64")
    
    incident["prevision"] = incident["prevision"].astype(str).apply(classify_prevision)
    
    return incident


# =============================================================================
# TASAS DE INCIDENCIA
# =============================================================================

def incidence_rates_and_counts(incident: pd.DataFrame, pop_long: pd.DataFrame) -> dict:
    """
    Calcula conteos y tasas de incidencia (2011-2023).
    
    Args:
        incident: DataFrame de casos incidentes
        pop_long: DataFrame de población INE
        
    Returns:
        Dict con:
        - cases_per_year: Serie de conteos por año
        - crude_incidence_per_year: Serie de tasas crudas
        - crude_incidence_per_year_by_sex: Dict {sexo: Serie}
        - age_specific_incidence: Serie MultiIndex (year, edad_grupo)
    """
    res = {}
    
    # Filtrar años de análisis
    pop_use = pop_long.loc[
        (pop_long['year'] >= INC_START_YEAR) & (pop_long['year'] <= MAX_YEAR)
    ].copy()
    
    pop_total = pop_use.groupby("year")["poblacion"].sum()
    pop_by_sex = pop_use.groupby(["year", "sexo"])['poblacion'].sum()

    inc_use = incident[
        (incident['year'] >= INC_START_YEAR) & (incident['year'] <= MAX_YEAR)
    ].copy()

    # Conteos por año
    cases_year = inc_use.groupby('year').size().rename('cases')
    res['cases_per_year'] = cases_year

    # Tasa cruda nacional
    rate_crude = (cases_year / pop_total * 1e5).dropna()
    res['crude_incidence_per_year'] = rate_crude

    # Tasa cruda por sexo
    if "sexo" in inc_use.columns:
        cases_year_sex = inc_use.groupby(['year', 'sexo']).size().rename('cases')
        rate_sex = (cases_year_sex / pop_by_sex * 1e5).dropna()
        res['crude_incidence_per_year_by_sex'] = {
            k: rate_sex.xs(k, level=1) 
            for k in rate_sex.index.get_level_values(1).unique()
        }

    # Tasa específica por edad
    pop_use.loc[:, 'edad_grupo'] = pop_use['Edad'].apply(to_age_group_5y)
    pop_age = pop_use.groupby(['year', 'edad_grupo'])['poblacion'].sum()
    cases_age = inc_use.groupby(['year', 'edad_grupo']).size().rename('cases')
    rate_age = (cases_age / pop_age * 1e5).dropna()
    res['age_specific_incidence'] = rate_age

    return res


# =============================================================================
# SERIES TEMPORALES CON INTERVALOS DE CONFIANZA
# =============================================================================

def incidence_age_sex_over_time(incident: pd.DataFrame, pop_long: pd.DataFrame) -> pd.DataFrame:
    """
    Tasas de incidencia por edad y sexo a lo largo del tiempo.
    
    Args:
        incident: DataFrame de casos incidentes
        pop_long: DataFrame de población
        
    Returns:
        DataFrame con columnas: year, edad_grupo, sexo, rate
    """
    inc = incident[
        (incident["year"] >= INC_START_YEAR) & (incident["year"] <= MAX_YEAR)
    ].copy()
    
    if inc.empty:
        return pd.DataFrame(columns=["year", "edad_grupo", "sexo", "rate"])

    pop_use = pop_long[
        (pop_long["year"] >= INC_START_YEAR) & (pop_long["year"] <= MAX_YEAR)
    ].copy()
    pop_use.loc[:, "edad_grupo"] = pop_use["Edad"].apply(lambda x: to_age_group_5y(x, TOP_OPEN))

    c = inc.groupby(["year", "edad_grupo", "sexo"]).size().rename("cases")
    p = pop_use.groupby(["year", "edad_grupo", "sexo"])["poblacion"].sum()
    rate = (c / p * 1e5).replace([np.inf, -np.inf], np.nan).dropna()
    
    df = rate.reset_index().rename(columns={0: "rate"})
    df["edad_grupo"] = pd.Categorical(
        df["edad_grupo"],
        categories=sort_age_index(df["edad_grupo"].unique()),
        ordered=True
    )
    
    return df.sort_values(["edad_grupo", "year", "sexo"])


def incidence_crude_timeseries_with_ci(incident: pd.DataFrame, pop_long: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """
    Series temporales de incidencia cruda con IC 95% por sexo y total.
    
    Args:
        incident: DataFrame de casos incidentes
        pop_long: DataFrame de población
        
    Returns:
        Dict {'Total', 'Female (F)', 'Male (M)'}: DataFrames con rate, lcl, ucl
    """
    years = list(range(INC_START_YEAR, MAX_YEAR + 1))
    
    # Total
    d_tot = incident.groupby("year").size().to_dict()
    p_tot = pop_long.groupby("year")["poblacion"].sum().to_dict()
    tot = {y: poisson_rate_ci(d_tot.get(y, 0), p_tot.get(y, np.nan)) for y in years}
    
    out = {"Total": _to_df_rate_ci(tot)}
    
    # Por sexo
    for sx, label in [("F", "Female (F)"), ("M", "Male (M)")]:
        d = incident[incident["sexo"] == sx].groupby("year").size().to_dict()
        p = pop_long.groupby(["year", "sexo"])["poblacion"].sum()
        p = {y: float(p.get((y, sx), np.nan)) for y in years}
        series = {y: poisson_rate_ci(d.get(y, 0), p.get(y, np.nan)) for y in years}
        out[label] = _to_df_rate_ci(series)
    
    return out


def incidence_asr_timeseries_with_ci(incident: pd.DataFrame, pop_long: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """
    Series temporales de incidencia ASR (WHO 80+) con IC aproximado.
    
    Args:
        incident: DataFrame de casos incidentes
        pop_long: DataFrame de población
        
    Returns:
        Dict {'Total', 'Female (F)', 'Male (M)'}: DataFrames con rate, lcl, ucl
    """
    years = list(range(INC_START_YEAR, MAX_YEAR + 1))
    who_w = who_world_standard_weights_collapsed()
    
    pop = pop_long.copy()
    pop["edad_grupo"] = pop["Edad"].apply(lambda x: to_age_group_5y(x, TOP_OPEN))

    out = {}
    
    # Total
    total_dict = {}
    for y in years:
        inc_y = incident[incident["year"] == y]
        if inc_y.empty:
            total_dict[y] = (np.nan, np.nan, np.nan)
            continue
        d_age = inc_y.groupby("edad_grupo").size()
        N_age = pop[pop["year"] == y].groupby("edad_grupo")["poblacion"].sum()
        asr, l, u = asr_and_ci_from_counts(d_age, N_age, who_w)
        total_dict[y] = (asr, l, u)
    out["Total"] = _to_df_rate_ci(total_dict)
    
    # Por sexo
    for sx, label in [("F", "Female (F)"), ("M", "Male (M)")]:
        dct = {}
        for y in years:
            inc_y = incident[(incident["year"] == y) & (incident["sexo"] == sx)]
            if inc_y.empty:
                dct[y] = (np.nan, np.nan, np.nan)
                continue
            d_age = inc_y.groupby("edad_grupo").size()
            N_age = pop[(pop["year"] == y) & (pop["sexo"] == sx)].groupby("edad_grupo")["poblacion"].sum()
            asr, l, u = asr_and_ci_from_counts(d_age, N_age, who_w)
            dct[y] = (asr, l, u)
        out[label] = _to_df_rate_ci(dct)
    
    return out


# =============================================================================
# UTILIDADES
# =============================================================================

def _to_df_rate_ci(dct: dict[int, tuple[float, float, float]]) -> pd.DataFrame:
    """Convierte dict de tuplas (rate, lcl, ucl) a DataFrame."""
    if not dct:
        return pd.DataFrame(columns=["rate", "lcl", "ucl"])
    df = pd.DataFrame.from_dict(dct, orient="index", columns=["rate", "lcl", "ucl"])
    df.index.name = "year"
    return df

