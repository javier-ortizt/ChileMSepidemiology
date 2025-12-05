"""
Análisis Estadístico
====================

Funciones para:
- Intervalos de confianza (Poisson, Pearson)
- Cálculo de ASR (Age-Standardized Rates) con IC
- Regresión Poisson para tendencias
- Correlación de Pearson con banda de confianza
"""

from __future__ import annotations
from math import atanh, tanh

import numpy as np
import pandas as pd
from scipy.stats import chi2, pearsonr, norm
import scipy.stats as stats
import statsmodels.api as sm
import statsmodels.formula.api as smf

from .who_standards import who_world_standard_weights_collapsed


# =============================================================================
# INTERVALOS DE CONFIANZA - POISSON
# =============================================================================

def poisson_count_ci(count: int, alpha: float = 0.05) -> tuple[float, float]:
    """
    Intervalo de confianza exacto para un conteo Poisson usando límites chi-cuadrado.
    
    Args:
        count: Número de casos observados
        alpha: Nivel de significancia (default 0.05 para IC 95%)
        
    Returns:
        Tupla (límite_inferior, límite_superior)
    """
    if count < 0:
        return (np.nan, np.nan)
    
    if count == 0:
        lower = 0.0
    else:
        lower = 0.5 * chi2.ppf(alpha / 2, 2 * count)
    
    upper = 0.5 * chi2.ppf(1 - alpha / 2, 2 * (count + 1))
    
    return (lower, upper)


def poisson_rate_ci(count: int, population: float, 
                    multiplier: float = 1e5, alpha: float = 0.05) -> tuple[float, float, float]:
    """
    Tasa por `multiplier` habitantes con IC exacto Poisson en el numerador.
    
    Args:
        count: Número de casos
        population: Población de referencia
        multiplier: Factor de multiplicación (default 100,000)
        alpha: Nivel de significancia
        
    Returns:
        Tupla (tasa, límite_inferior, límite_superior)
    """
    if population <= 0 or pd.isna(population):
        return (np.nan, np.nan, np.nan)
    
    rate = count / population * multiplier
    lo_c, hi_c = poisson_count_ci(int(count), alpha)
    lcl = lo_c / population * multiplier
    ucl = hi_c / population * multiplier
    
    return (rate, lcl, ucl)


# =============================================================================
# INTERVALO DE CONFIANZA - PEARSON
# =============================================================================

def pearson_ci(r: float, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """
    IC para correlación de Pearson usando transformación Z de Fisher.
    
    Args:
        r: Coeficiente de correlación de Pearson
        n: Tamaño de muestra
        alpha: Nivel de significancia
        
    Returns:
        Tupla (límite_inferior, límite_superior)
    """
    if n < 4:
        return (np.nan, np.nan)
    
    z = atanh(r)
    se = 1 / np.sqrt(n - 3)
    z_crit = norm.ppf(1 - alpha / 2)
    lo, hi = z - z_crit * se, z + z_crit * se
    
    return (tanh(lo), tanh(hi))


# =============================================================================
# CÁLCULO DE ASR (AGE-STANDARDIZED RATES)
# =============================================================================

def asr_from_rates(rate_age: pd.Series, weights: dict) -> pd.Series:
    """
    Calcula ASR WHO por año desde tasas específicas por edad.
    
    Re-normaliza pesos según grupos de edad disponibles.
    
    Args:
        rate_age: Serie con MultiIndex (year, edad_grupo) y tasas por 100,000
        weights: Dict de pesos WHO por grupo de edad
        
    Returns:
        Serie indexada por año con ASR
    """
    if rate_age is None or rate_age.empty:
        return pd.Series(dtype=float, name="WHO_ASR_per_100k")
    
    df = rate_age.reset_index(name="rate")
    df = df[df['edad_grupo'].isin(weights.keys())].copy()
    
    if df.empty:
        return pd.Series(dtype=float, name="WHO_ASR_per_100k")
    
    df['w'] = df['edad_grupo'].map(weights)

    def _asr(g):
        w = g['w'].to_numpy()
        r = g['rate'].to_numpy()
        w = w / w.sum() if w.sum() else w
        return float((r * w).sum())

    asr = df.groupby('year', sort=True)[['rate', 'w']].apply(lambda g: _asr(g))
    asr.name = 'WHO_ASR_per_100k'
    
    return asr


def asr_from_rates_by_group(rate_age_df: pd.DataFrame, weights: dict, group_col: str) -> pd.Series:
    """
    Calcula ASR por año y grupo adicional (ej: región, sexo).
    
    Args:
        rate_age_df: DataFrame con columnas year, edad_grupo, group_col, rate
        weights: Dict de pesos WHO
        group_col: Nombre de la columna de agrupación adicional
        
    Returns:
        Serie con MultiIndex (year, group_col) y ASR
    """
    if rate_age_df.empty:
        return pd.Series(dtype=float, name="WHO_ASR_per_100k")

    def _one(g):
        g = g[g['edad_grupo'].isin(weights.keys())]
        if g.empty:
            return np.nan
        w = g['edad_grupo'].map(weights).to_numpy()
        r = g['rate'].to_numpy()
        w = w / w.sum() if w.sum() else w
        return float((r * w).sum())

    out = (
        rate_age_df
        .groupby(['year', group_col], sort=True)
        .apply(_one, include_groups=False)
        .dropna()
    )
    out.name = "WHO_ASR_per_100k"
    
    return out


def asr_by_sex_from_age_rates(rate_age_sex_df: pd.DataFrame, weights: dict) -> dict[str, pd.Series]:
    """
    Calcula ASR por año y sexo.
    
    Args:
        rate_age_sex_df: DataFrame con columnas year, edad_grupo, sexo, rate
        weights: Dict de pesos WHO
        
    Returns:
        Dict {'F': Series(year), 'M': Series(year)} con ASR
    """
    if rate_age_sex_df is None or rate_age_sex_df.empty:
        return {"F": pd.Series(dtype=float), "M": pd.Series(dtype=float)}
    
    out = {}
    for sx in ["F", "M"]:
        sub = rate_age_sex_df[rate_age_sex_df["sexo"] == sx][["year", "edad_grupo", "rate"]]
        if sub.empty:
            out[sx] = pd.Series(dtype=float, name="WHO_ASR_per_100k")
            continue
        out[sx] = asr_from_rates(sub.set_index(["year", "edad_grupo"])["rate"], weights)
    
    return out


def asr_and_ci_from_counts(cases_by_age: pd.Series, pop_by_age: pd.Series,
                           weights: dict, multiplier: float = 1e5, 
                           alpha: float = 0.05) -> tuple[float, float, float]:
    """
    Calcula ASR con IC aproximado desde conteos y población por edad.
    
    Fórmula:
        ASR = sum_i w_i * (d_i / N_i) * multiplier
        Var(ASR) ≈ sum_i w_i^2 * (d_i / N_i^2) * multiplier^2
        IC = ASR ± z * sqrt(Var)
    
    Args:
        cases_by_age: Serie de casos por grupo de edad
        pop_by_age: Serie de población por grupo de edad
        weights: Dict de pesos WHO
        multiplier: Factor de multiplicación (default 100,000)
        alpha: Nivel de significancia
        
    Returns:
        Tupla (ASR, límite_inferior, límite_superior)
    """
    df = pd.DataFrame({"d": cases_by_age, "N": pop_by_age}).dropna()
    df = df[df.index.isin(weights.keys())]
    
    if df.empty:
        return (np.nan, np.nan, np.nan)
    
    # Pesos re-normalizados
    w = pd.Series({k: weights[k] for k in df.index})
    w = w / w.sum() if w.sum() else w
    
    # Tasas por edad
    rate_i = (df["d"] / df["N"]).fillna(0.0)
    
    # ASR
    asr = float((w * rate_i).sum() * multiplier)
    
    # Varianza aproximada
    var = float(((w**2) * (df["d"] / (df["N"]**2))).sum() * (multiplier**2))
    se = np.sqrt(var)
    
    z = norm.ppf(1 - alpha / 2)
    
    return (asr, asr - z * se, asr + z * se)


# =============================================================================
# REGRESIÓN POISSON PARA TENDENCIAS
# =============================================================================

def _fit_poisson(formula: str, data: pd.DataFrame, offset_col: str):
    """Ajusta modelo Poisson GLM con offset logarítmico."""
    return smf.glm(
        formula=formula,
        data=data,
        family=sm.families.Poisson(),
        offset=np.log(data[offset_col]),
    ).fit(cov_type="HC0")


def _poisson_print(model, label: str) -> None:
    """Imprime resultados de modelo Poisson."""
    if 'year' not in model.params.index:
        print(f"\n{label}\n  [WARNING] El modelo no tiene coeficiente para 'year'.")
        return
    
    b = model.params['year']
    se = model.bse['year']
    irr = float(np.exp(b))
    lcl = float(np.exp(b - 1.96 * se))
    ucl = float(np.exp(b + 1.96 * se))
    p = float(model.pvalues['year']) if 'year' in model.pvalues.index else np.nan
    
    pearson = float((model.resid_pearson**2).sum())
    disp = pearson / model.df_resid if model.df_resid > 0 else np.nan
    
    print(f"\n{label}")
    print(f"  IRR per year = {irr:.4f} (95%CI {lcl:.4f} – {ucl:.4f}), p = {p:.3g}")
    print(f"  Dispersion (Pearson/df) = {disp:.3f}")


def poisson_trends_incidence(incident: pd.DataFrame, pop_long: pd.DataFrame,
                             year_min: int, year_max: int) -> None:
    """
    Analiza tendencia temporal de incidencia usando regresión Poisson.
    
    Args:
        incident: DataFrame de casos incidentes
        pop_long: DataFrame de población
        year_min, year_max: Rango de años
    """
    inc = incident[(incident['year'] >= year_min) & (incident['year'] <= year_max)].copy()
    
    if inc.empty:
        print("\n[WARNING] No hay incidentes en el rango para análisis Poisson.")
        return
    
    pop_use = pop_long[(pop_long['year'] >= year_min) & (pop_long['year'] <= year_max)].copy()

    # Total (crudo)
    c_tot = inc.groupby('year').size().rename('count').reset_index()
    p_tot = pop_use.groupby('year')['poblacion'].sum().rename('pop').reset_index()
    d_tot = c_tot.merge(p_tot, on='year', how='inner')
    d_tot = d_tot[d_tot['pop'] > 0]
    
    m_tot = _fit_poisson('count ~ year', d_tot, 'pop')
    _poisson_print(m_tot, f"POISSON (Incidencia TOTAL {year_min}–{year_max}) — tendencia cruda")


def poisson_trends_prevalence(prev_obj: dict, pop_long: pd.DataFrame,
                              year_min: int, year_max: int) -> None:
    """
    Analiza tendencia temporal de prevalencia usando regresión Poisson.
    
    Args:
        prev_obj: Dict con conteos de prevalencia
        pop_long: DataFrame de población
        year_min, year_max: Rango de años
    """
    pop_use = pop_long[(pop_long['year'] >= year_min) & (pop_long['year'] <= year_max)].copy()
    prev_tot = prev_obj.get('prev_counts', pd.Series(dtype=int))
    prev_tot = prev_tot.loc[lambda s: (s.index >= year_min) & (s.index <= year_max)]
    
    if prev_tot.empty:
        print("\n[WARNING] No hay prevalentes para análisis Poisson.")
        return
    
    d_tot = prev_tot.rename('count').reset_index().rename(columns={'index': 'year'})
    p_tot = pop_use.groupby('year')['poblacion'].sum().rename('pop').reset_index()
    d_tot = d_tot.merge(p_tot, on='year', how='inner')
    d_tot = d_tot[d_tot['pop'] > 0]
    
    m_tot = _fit_poisson('count ~ year', d_tot, 'pop')
    _poisson_print(m_tot, f"POISSON (Prevalencia Vivos TOTAL {year_min}–{year_max}) — tendencia cruda")


# =============================================================================
# BANDA DE CONFIANZA PARA REGRESIÓN OLS (PARA GRÁFICOS)
# =============================================================================

def ols_confidence_band(x: np.ndarray, y: np.ndarray, 
                        xfit: np.ndarray, alpha: float = 0.05):
    """
    Calcula banda de confianza para predicción media de OLS y ~ x.
    
    Args:
        x, y: Datos observados
        xfit: Valores de x para los que calcular la predicción
        alpha: Nivel de significancia
        
    Returns:
        Tupla (yfit, lcl, ucl) - predicción y límites de confianza
    """
    x1 = np.column_stack([np.ones_like(x), x])
    beta, *_ = np.linalg.lstsq(x1, y, rcond=None)
    yfit = beta[0] + beta[1] * xfit
    yhat = beta[0] + beta[1] * x
    
    n = len(x)
    sse = np.sum((y - yhat)**2)
    se2 = sse / (n - 2)
    xbar = np.mean(x)
    sxx = np.sum((x - xbar)**2)
    
    tcrit = stats.t.ppf(1 - alpha / 2, df=n - 2)
    se_mean = np.sqrt(se2 * (1/n + (xfit - xbar)**2 / sxx))
    
    lcl = yfit - tcrit * se_mean
    ucl = yfit + tcrit * se_mean
    
    return yfit, lcl, ucl

