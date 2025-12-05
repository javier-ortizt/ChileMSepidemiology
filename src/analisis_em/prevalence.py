"""
Cálculo de Prevalencia
======================

Funciones para el cálculo de prevalencia de EM:
- Base de prevalentes (beneficiarios únicos vivos)
- Conteos y tasas por año, sexo, edad, región
- Series temporales con intervalos de confianza
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .config import MIN_YEAR, MAX_YEAR, TOP_OPEN
from .utils import (
    find_column, classify_prevision, normalize_sex,
    to_age_group_5y, to_age_group_10y, sort_age_index, age_group_sortkey
)
from .statistics import (
    asr_from_rates, asr_from_rates_by_group, poisson_rate_ci, asr_and_ci_from_counts
)
from .who_standards import who_world_standard_weights_collapsed, who_world_standard_weights_10y
from .incidence import extract_year_column


# =============================================================================
# BASE DE PREVALENCIA
# =============================================================================

def build_prevalence_base(df_cases: pd.DataFrame) -> pd.DataFrame:
    """
    Construye base de prevalencia con un registro único por beneficiario.
    
    Un caso es considerado VIVO (prevalente) si death_is_null == True
    (es decir, no tiene fecha de fallecimiento registrada).
    
    Args:
        df_cases: DataFrame con todos los registros de casos
        
    Returns:
        DataFrame con beneficiarios únicos:
        - diag_year: Año de diagnóstico
        - death_is_null: True si está vivo (sin fecha de fallecimiento)
        - sexo: M/F
        - edad: Edad al diagnóstico
        - edad_grupo: Grupo quinquenal
        - region_code: Código de región
        - prevision: FONASA/ISAPRE/FFAA
        - prevision_raw: Valor original de previsión
    """
    df = df_cases.copy()

    # Año de diagnóstico
    year_series = extract_year_column(df)
    if year_series is not None:
        df["diag_year"] = year_series
    else:
        f_diag = find_column(df, ["fecha confirmacion diagnostica", "fecha confirmación diagnóstica"])
        if not f_diag:
            f_diag = find_column(df, ["fecha creacion caso ges", "fecha creación caso ges"])
        df["diag_year"] = (
            pd.to_datetime(df[f_diag], errors="coerce", dayfirst=True).dt.year 
            if f_diag else np.nan
        )

    # Fecha de fallecimiento -> indicador de vivo
    f_dead = find_column(df, ["fecha fallecimiento beneficiario", "fecha fallecimiento"])
    death_dt = (
        pd.to_datetime(df[f_dead], errors="coerce", dayfirst=True) 
        if f_dead else pd.Series(pd.NaT, index=df.index)
    )
    df["death_is_null"] = death_dt.isna()

    # Identificar columnas relevantes
    id_col = find_column(df, ["cod beneficiario", "id beneficiario", "id", "codigo beneficiario", "rut", "n°", "n "])
    sexo_col = find_column(df, ["sexo beneficiario", "sexo"])
    edad_col = find_column(df, ["edad (anos) creacion caso ges", "edad (años) creacion caso ges", "edad"])
    region_code_col = find_column(df, ["cod. deis region beneficiario", "cod deis region beneficiario", "cod deis region"])
    prev_col = find_column(df, ["prevision"])

    # Seleccionar columnas
    keep = [c for c in [id_col, "diag_year", "death_is_null", sexo_col, edad_col, region_code_col, prev_col] if c]
    base = df[keep].dropna(subset=[id_col, "diag_year"]).copy()
    
    # Mantener primer registro por beneficiario
    base.sort_values([id_col, "diag_year"], inplace=True)
    base = base.drop_duplicates(subset=[id_col], keep="first")
    
    # Renombrar columnas
    if sexo_col:
        base.rename(columns={sexo_col: "sexo"}, inplace=True)
    if edad_col:
        base.rename(columns={edad_col: "edad"}, inplace=True)
    if region_code_col:
        base.rename(columns={region_code_col: "region_code"}, inplace=True)
        base["region_code"] = pd.to_numeric(base["region_code"], errors="coerce").astype("Int64")
    if prev_col:
        base.rename(columns={prev_col: "prevision_raw"}, inplace=True)
        base["prevision_raw"] = base["prevision_raw"].astype(str)
        base["prevision"] = base["prevision_raw"].apply(classify_prevision)

    # Normalizar valores
    base["sexo"] = base["sexo"].apply(normalize_sex)
    base["edad_grupo"] = base["edad"].apply(lambda x: to_age_group_5y(x, TOP_OPEN))
    
    return base


# =============================================================================
# CONTEOS DE PREVALENCIA
# =============================================================================

def prevalence_counts(base: pd.DataFrame, pop_long: pd.DataFrame) -> dict:
    """
    Construye conteos de prevalentes por año, sexo, edad, región y previsión.
    
    Un caso es prevalente en año Y si:
    - Está vivo (death_is_null == True)
    - Fue diagnosticado hasta Y (diag_year <= Y)
    
    Args:
        base: DataFrame de base de prevalencia
        pop_long: DataFrame de población INE
        
    Returns:
        Dict con:
        - years: Lista de años
        - prev_counts: Serie de conteos totales por año
        - prev_counts_sex: Dict {año: Serie por sexo}
        - prev_counts_age: Dict {año: Serie por edad}
        - prev_counts_region: Dict {año: Serie por región}
        - prev_counts_prev: Dict {año: Serie por previsión}
    """
    years_pop = sorted(pop_long["year"].unique())
    years = [y for y in years_pop if MIN_YEAR <= y <= MAX_YEAR]
    base_alive = base[base["death_is_null"]].copy()

    prev_counts = {}
    prev_counts_sex = {}
    prev_counts_age = {}
    prev_counts_region = {}
    prev_counts_prev = {}

    for y in years:
        prev_y = base_alive.loc[base_alive["diag_year"] <= y]
        prev_counts[y] = len(prev_y)
        
        if "sexo" in prev_y.columns:
            prev_counts_sex[y] = prev_y.groupby("sexo").size()
        if "edad_grupo" in prev_y.columns:
            prev_counts_age[y] = prev_y.groupby("edad_grupo").size()
        if "region_code" in prev_y.columns:
            prev_counts_region[y] = prev_y.groupby("region_code").size()
        if "prevision" in prev_y.columns:
            prev_counts_prev[y] = prev_y.groupby("prevision").size()

    return {
        "years": years,
        "prev_counts": pd.Series(prev_counts, name="cases"),
        "prev_counts_sex": prev_counts_sex,
        "prev_counts_age": prev_counts_age,
        "prev_counts_region": prev_counts_region,
        "prev_counts_prev": prev_counts_prev,
    }


# =============================================================================
# TASAS DE PREVALENCIA
# =============================================================================

def prevalence_rates(prev_obj: dict, pop_long: pd.DataFrame) -> dict:
    """
    Calcula tasas de prevalencia desde conteos.
    
    Args:
        prev_obj: Dict de conteos de prevalencia_counts()
        pop_long: DataFrame de población INE
        
    Returns:
        Dict con:
        - crude_prevalence_per_year: Serie de tasas crudas
        - crude_prevalence_per_year_by_sex: Dict {sexo: Serie}
        - age_specific_prevalence: Serie MultiIndex (year, edad_grupo)
        - crude_prevalence_per_year_by_region: Serie MultiIndex (year, region_code)
        - who_asr_per_year: Serie ASR por año
    """
    res = {}
    pop_use = pop_long.loc[
        (pop_long['year'] >= MIN_YEAR) & (pop_long['year'] <= MAX_YEAR)
    ].copy()

    # Tasa cruda nacional
    pop_total = pop_use.groupby("year")["poblacion"].sum()
    prev_crude = (prev_obj["prev_counts"] / pop_total * 1e5).dropna()
    res["crude_prevalence_per_year"] = prev_crude

    # Por sexo (cruda)
    if prev_obj["prev_counts_sex"]:
        pop_by_sex = pop_use.groupby(["year", "sexo"])['poblacion'].sum()
        records = []
        for y, s in prev_obj["prev_counts_sex"].items():
            for sex, n in s.items():
                records.append((y, str(sex), int(n)))
        df_sex = pd.DataFrame(records, columns=["year", "sexo", "cases"])
        rate_sex = (df_sex.set_index(["year", "sexo"])["cases"] / pop_by_sex * 1e5).dropna()
        res["crude_prevalence_per_year_by_sex"] = {
            k: rate_sex.xs(k, level=1) 
            for k in rate_sex.index.get_level_values(1).unique()
        }

    # Específica por edad (nacional)
    pop_use.loc[:, "edad_grupo"] = pop_use["Edad"].apply(lambda x: to_age_group_5y(x, TOP_OPEN))
    pop_age = pop_use.groupby(["year", "edad_grupo"])['poblacion'].sum()
    
    records = []
    for y, s in prev_obj["prev_counts_age"].items():
        for ag, n in s.items():
            records.append((y, str(ag), int(n)))
    
    if records:
        df_age = pd.DataFrame(records, columns=["year", "edad_grupo", "cases"])
        df_age = df_age.set_index(["year", "edad_grupo"]).sort_index()
        rate_age = (df_age["cases"] / pop_age * 1e5).dropna()
        res["age_specific_prevalence"] = rate_age
        
        # WHO ASR (nacional)
        res["who_asr_per_year"] = asr_from_rates(rate_age, who_world_standard_weights_collapsed())
    else:
        res["age_specific_prevalence"] = pd.Series(dtype=float)
        res["who_asr_per_year"] = pd.Series(dtype=float)

    # Por región (cruda)
    if prev_obj["prev_counts_region"]:
        pop_region = pop_use.groupby(["year", "region_code"])['poblacion'].sum()
        records = []
        for y, s in prev_obj["prev_counts_region"].items():
            for regcode, n in s.items():
                records.append((y, int(regcode) if pd.notna(regcode) else np.nan, int(n)))
        df_reg = pd.DataFrame(records, columns=["year", "region_code", "cases"])
        df_reg = df_reg.dropna(subset=["region_code"])
        df_reg["region_code"] = df_reg["region_code"].astype(int)
        rate_reg = (df_reg.set_index(["year", "region_code"])["cases"] / pop_region).mul(1e5).dropna()
        res["crude_prevalence_per_year_by_region"] = rate_reg

    return res


# =============================================================================
# PREVALENCIA POR EDAD Y REGIÓN (PARA ASR REGIONAL)
# =============================================================================

def prevalence_by_region_age_for_asr(base_prev: pd.DataFrame, pop_long: pd.DataFrame) -> pd.DataFrame:
    """
    Construye tasas de prevalencia por edad y región para cálculo de ASR regional.
    
    Args:
        base_prev: DataFrame de base de prevalencia
        pop_long: DataFrame de población INE
        
    Returns:
        DataFrame con columnas: year, region_code, edad_grupo, rate
    """
    if base_prev is None or base_prev.empty:
        return pd.DataFrame(columns=["year", "region_code", "edad_grupo", "rate"])

    alive = base_prev[base_prev["death_is_null"]].copy()
    alive = alive.dropna(subset=["diag_year", "region_code", "edad_grupo"])
    alive["region_code"] = alive["region_code"].astype(int)

    years = [y for y in sorted(pop_long["year"].unique()) if MIN_YEAR <= y <= MAX_YEAR]
    out = []

    # Pre-calcular población por año × región × edad
    pop = pop_long[pop_long["year"].isin(years)].copy()
    pop["edad_grupo"] = pop["Edad"].apply(lambda x: to_age_group_5y(x, TOP_OPEN))
    pop_grp = pop.groupby(["year", "region_code", "edad_grupo"])["poblacion"].sum()

    for y in years:
        sub = alive[alive["diag_year"] <= y]
        if sub.empty:
            continue
        cases = sub.groupby(["region_code", "edad_grupo"]).size().rename("cases")
        try:
            p = pop_grp.xs(y, level=0)
        except KeyError:
            continue
        common = cases.index.intersection(p.index)
        if len(common) == 0:
            continue
        rates = (cases.loc[common] / p.loc[common] * 1e5).replace([np.inf, -np.inf], np.nan).dropna()
        if not rates.empty:
            tmp = rates.reset_index()
            tmp.insert(0, "year", y)
            tmp.rename(columns={0: "rate"}, inplace=True)
            tmp.columns = ["year", "region_code", "edad_grupo", "rate"]
            out.append(tmp)

    if not out:
        return pd.DataFrame(columns=["year", "region_code", "edad_grupo", "rate"])
    
    df_rates = pd.concat(out, ignore_index=True)
    df_rates["edad_grupo"] = pd.Categorical(
        df_rates["edad_grupo"],
        categories=sorted(df_rates["edad_grupo"].unique(), key=age_group_sortkey),
        ordered=True
    )
    return df_rates


# =============================================================================
# SERIES TEMPORALES CON INTERVALOS DE CONFIANZA
# =============================================================================

def prevalence_crude_timeseries_with_ci(base_prev: pd.DataFrame, pop_long: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """
    Series temporales de prevalencia cruda con IC 95%.
    
    Args:
        base_prev: DataFrame de base de prevalencia
        pop_long: DataFrame de población
        
    Returns:
        Dict {'Total', 'Female (F)', 'Male (M)'}: DataFrames con rate, lcl, ucl
    """
    years = list(range(MIN_YEAR, MAX_YEAR + 1))
    alive = base_prev[base_prev["death_is_null"]].copy()

    out = {}
    
    # Total
    tot = {}
    for y in years:
        d = int((alive["diag_year"] <= y).sum())
        N = float(pop_long[pop_long["year"] == y]["poblacion"].sum())
        tot[y] = poisson_rate_ci(d, N, 1e5)
    out["Total"] = _to_df_rate_ci(tot)
    
    # Por sexo
    for sx, label in [("F", "Female (F)"), ("M", "Male (M)")]:
        dct = {}
        alive_sx = alive[alive["sexo"] == sx]
        for y in years:
            d = int((alive_sx["diag_year"] <= y).sum())
            N = float(pop_long[(pop_long["year"] == y) & (pop_long["sexo"] == sx)]["poblacion"].sum())
            dct[y] = poisson_rate_ci(d, N, 1e5)
        out[label] = _to_df_rate_ci(dct)
    
    return out


def prevalence_asr_timeseries_with_ci(base_prev: pd.DataFrame, pop_long: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """
    Series temporales de prevalencia ASR (WHO 80+) con IC aproximado.
    
    Args:
        base_prev: DataFrame de base de prevalencia
        pop_long: DataFrame de población
        
    Returns:
        Dict {'Total', 'Female (F)', 'Male (M)'}: DataFrames con rate, lcl, ucl
    """
    who_w = who_world_standard_weights_collapsed()
    years = list(range(MIN_YEAR, MAX_YEAR + 1))
    alive = base_prev[base_prev["death_is_null"]].copy()
    
    pop = pop_long.copy()
    pop["edad_grupo"] = pop["Edad"].apply(lambda x: to_age_group_5y(x, TOP_OPEN))

    out = {}
    
    # Total
    tot = {}
    for y in years:
        sub = alive[alive["diag_year"] <= y]
        d_age = sub.groupby("edad_grupo").size()
        N_age = pop[pop["year"] == y].groupby("edad_grupo")["poblacion"].sum()
        asr, l, u = asr_and_ci_from_counts(d_age, N_age, who_w)
        tot[y] = (asr, l, u)
    out["Total"] = _to_df_rate_ci(tot)
    
    # Por sexo
    for sx, label in [("F", "Female (F)"), ("M", "Male (M)")]:
        dct = {}
        alive_sx = alive[alive["sexo"] == sx]
        for y in years:
            sub = alive_sx[alive_sx["diag_year"] <= y]
            d_age = sub.groupby("edad_grupo").size()
            N_age = pop[(pop["year"] == y) & (pop["sexo"] == sx)].groupby("edad_grupo")["poblacion"].sum()
            asr, l, u = asr_and_ci_from_counts(d_age, N_age, who_w)
            dct[y] = (asr, l, u)
        out[label] = _to_df_rate_ci(dct)
    
    return out


def prevalence_age_sex_over_time(base_prev: pd.DataFrame, pop_long: pd.DataFrame) -> pd.DataFrame:
    """
    Tasas de prevalencia por edad y sexo a lo largo del tiempo.
    
    Args:
        base_prev: DataFrame de base de prevalencia
        pop_long: DataFrame de población
        
    Returns:
        DataFrame con columnas: year, edad_grupo, sexo, rate
    """
    alive = base_prev[base_prev["death_is_null"]].copy()
    weights_years = range(MIN_YEAR, MAX_YEAR + 1)
    
    pop = pop_long[pop_long["year"].isin(weights_years)].copy()
    pop["edad_grupo"] = pop["Edad"].apply(lambda x: to_age_group_5y(x, TOP_OPEN))
    p = pop.groupby(["year", "edad_grupo", "sexo"])["poblacion"].sum()

    records = []
    for y in weights_years:
        sub = alive[alive["diag_year"] <= y]
        if sub.empty:
            continue
        c2 = sub.groupby(["edad_grupo", "sexo"]).size().rename("cases").reset_index()
        c2["year"] = y
        c2 = c2.set_index(["year", "edad_grupo", "sexo"])["cases"]
        aligned = (c2 / p).mul(1e5).replace([np.inf, -np.inf], np.nan).dropna()
        if not aligned.empty:
            dfy = aligned.reset_index().rename(columns={0: "rate"})
            dfy.columns = ["year", "edad_grupo", "sexo", "rate"]
            records.append(dfy)
    
    if not records:
        return pd.DataFrame(columns=["year", "edad_grupo", "sexo", "rate"])
    
    df = pd.concat(records, ignore_index=True)
    df["edad_grupo"] = pd.Categorical(
        df["edad_grupo"],
        categories=sort_age_index(df["edad_grupo"].unique()),
        ordered=True
    )
    return df.sort_values(["year", "edad_grupo", "sexo"])


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

