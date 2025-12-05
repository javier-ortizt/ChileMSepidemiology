"""
Función Principal
=================

Orquesta el análisis completo de EM:
1. Carga de datos
2. Cálculo de incidencia y prevalencia
3. Análisis por sistema previsional
4. Generación de gráficos y mapas
5. Exportación de resultados
"""

from __future__ import annotations
from pathlib import Path
import warnings

import pandas as pd
import numpy as np

from .config import (
    EXCEL_FILE, INE_POP_FILE, PREV_REGION_FILE, BEN_SYS_FILE, SHAPEFILE_PATH,
    MIN_YEAR, MAX_YEAR, INC_START_YEAR,
    get_output_dirs
)
from .utils import find_column, normalize_sex, to_age_group_5y, age_group_sortkey, TOP_OPEN
from .io_data import (
    read_bdem, load_ine_population,
    load_population_by_prevision_region, load_population_by_system_age
)
from .incidence import compute_incidence_from_bdem, incidence_rates_and_counts
from .prevalence import build_prevalence_base, prevalence_counts, prevalence_rates
from .statistics import asr_from_rates
from .who_standards import who_world_standard_weights_collapsed
from .insurance import (
    crude_prevalence_by_system, debug_crude_prevalence_by_system, who_asr_by_system,
    plot_isapre_detail_by_year
)
from .correlation import pearson_prevalence_vs_latitude
from .plotting import (
    save_line, save_bar, save_two_series_line,
    plot_incidence_by_year_crude_and_asr_with_ci,
    plot_prevalence_by_year_crude_and_asr_with_ci,
    plot_new_cases_by_year_by_system,
    plot_mean_age_at_incidence,
    plot_sex_ratio_and_percent_prevalent,
    plot_prevalence_by_system_over_time,
    prevalence_age_hist_with_totals
)
from .maps import plot_prevalence_map_by_region


# =============================================================================
# REPORTE DESCRIPTIVO
# =============================================================================

def print_overall_report(df: pd.DataFrame) -> None:
    """Imprime estadísticas descriptivas básicas del dataset."""
    print("\n=== REPORTE DESCRIPTIVO ===")
    
    col_id = find_column(df, ["cod beneficiario", "id beneficiario", "id", "codigo beneficiario", "rut", "n°", "n "])
    col_age = find_column(df, ["edad (anos) creacion caso ges", "edad (años) creacion caso ges", "edad"])
    col_sex = find_column(df, ["sexo beneficiario", "sexo"])

    total_rows = len(df)
    total_unique = df[col_id].nunique(dropna=True) if col_id in df.columns else total_rows
    
    print(f"Total filas: {total_rows}")
    print(f"Beneficiarios únicos: {total_unique}")

    # Edad
    if col_age and col_age in df.columns:
        ages = pd.to_numeric(df[col_age], errors="coerce")
        grp = ages.dropna().apply(lambda x: to_age_group_5y(x, TOP_OPEN))
        cnt = pd.Series(grp).value_counts().sort_index(key=lambda idx: [age_group_sortkey(i) for i in idx])
        pct = (cnt / cnt.sum() * 100).round(2)
        
        print("\n> Edad")
        print("Conteo por grupo de edad:")
        print(cnt.to_string())
        print("\nPorcentaje por grupo de edad (%):")
        print(pct.to_string())
        print(f"\nEdad media: {ages.mean():.2f}")

    # Sexo
    if col_sex and col_sex in df.columns:
        sx = df[col_sex].apply(normalize_sex)
        sxc = sx.value_counts(dropna=True)
        sxp = (sxc / sxc.sum() * 100).round(2)
        
        print("\n> Sexo")
        print("Conteo por sexo:")
        print(sxc.to_string())
        print("\nPorcentaje por sexo (%):")
        print(sxp.to_string())
    
    print("===========================\n")


# =============================================================================
# CÁLCULO Y REPORTE DE ASR DE INCIDENCIA
# =============================================================================

def compute_and_print_incidence_who_asr(inc_obj: dict, outdir: Path, gdir: Path) -> pd.Series:
    """Calcula y reporta ASR WHO de incidencia."""
    inc_age = inc_obj.get("age_specific_incidence", pd.Series(dtype=float))
    
    if inc_age is None or inc_age.empty:
        print("\n[WARNING] No se encontraron tasas específicas por edad para calcular WHO ASR.")
        return pd.Series(dtype=float)

    asr_inc = asr_from_rates(inc_age, who_world_standard_weights_collapsed())
    
    if asr_inc is None or asr_inc.empty:
        print("\n[WARNING] No se pudo calcular WHO ASR de incidencia.")
        return pd.Series(dtype=float)

    asr_inc = asr_inc.reindex(range(INC_START_YEAR, MAX_YEAR + 1))
    
    # Exportar
    asr_path = outdir / "incidence_who_asr_per_year.csv"
    asr_inc.round(2).to_csv(asr_path, encoding="utf-8-sig")
    
    save_line(asr_inc, "WHO-standardized incidence (80+) (per 100,000) – 2011–2023",
              gdir / "incidence_who_asr_per_year.png")
    
    print("\nIncidencia WHO-estandarizada (80+) (por 100,000) 2011–2023")
    print("[Nota] FFAA<2010 reasignados a 2010 (baseline); 2010 excluido del reporte.")
    print(asr_inc.round(2).to_string())
    
    return asr_inc


# =============================================================================
# FUNCIÓN PRINCIPAL
# =============================================================================

def main() -> None:
    """Ejecuta el análisis completo de EM."""
    
    # Crear directorios de salida
    outdir, gdir = get_output_dirs()

    # ========================================
    # 1. CARGA DE DATOS
    # ========================================
    print("=" * 60)
    print("CARGANDO DATOS...")
    print("=" * 60)
    
    df = read_bdem(EXCEL_FILE)
    pop_long = load_ine_population(INE_POP_FILE)
    sys_pop_age = load_population_by_system_age(BEN_SYS_FILE, year=MAX_YEAR)

    try:
        sys_pop_df = load_population_by_prevision_region(PREV_REGION_FILE, pop_long=pop_long, year=MAX_YEAR)
    except Exception as e:
        warnings.warn(f"[WARN] No se cargó población por previsión: {e}")
        sys_pop_df = None

    # ========================================
    # 2. REPORTE DESCRIPTIVO
    # ========================================
    print_overall_report(df)

    # ========================================
    # 3. ANÁLISIS DE INCIDENCIA
    # ========================================
    print("=" * 60)
    print("CALCULANDO INCIDENCIA...")
    print("=" * 60)
    
    incident = compute_incidence_from_bdem(df)
    inc = incidence_rates_and_counts(incident, pop_long)

    # Casos nuevos por año
    new_cases = inc["cases_per_year"].reindex(range(INC_START_YEAR, MAX_YEAR + 1)).fillna(0).astype(int)
    new_cases.to_csv(outdir / "new_cases_per_year.csv", encoding="utf-8-sig")
    save_bar(new_cases, "New cases per year (count) – 2011–2023", gdir / "new_cases_per_year.png")
    print("\nCASOS NUEVOS por año (conteo) 2011–2023:")
    print(new_cases.to_string())

    # Incidencia cruda
    inc_crude = inc["crude_incidence_per_year"].round(2).reindex(range(INC_START_YEAR, MAX_YEAR + 1))
    inc_crude.to_csv(outdir / "crude_incidence_per_year.csv", encoding="utf-8-sig")
    save_line(inc_crude, "Crude incidence per year (per 100,000) – 2011–2023", gdir / "crude_incidence_per_year.png")
    print("\nINCIDENCIA CRUDA (por 100,000) 2011–2023:")
    print(inc_crude.to_string())

    # Incidencia específica por edad y ASR
    inc_age = inc.get("age_specific_incidence", pd.Series(dtype=float))
    if inc_age is not None and not inc_age.empty:
        inc_age = inc_age[inc_age.index.get_level_values(0).isin(range(INC_START_YEAR, MAX_YEAR + 1))]
        inc_age.to_csv(outdir / "age_specific_incidence.csv", header=True, encoding="utf-8-sig")
        asr_inc = compute_and_print_incidence_who_asr(inc, outdir, gdir)
    else:
        asr_inc = pd.Series(dtype=float)

    # ========================================
    # 4. ANÁLISIS DE PREVALENCIA
    # ========================================
    print("\n" + "=" * 60)
    print("CALCULANDO PREVALENCIA...")
    print("=" * 60)
    
    base_prev = build_prevalence_base(df)
    
    # Prevalencia cruda por sistema
    crude = crude_prevalence_by_system(base_prev, sys_pop_age, MAX_YEAR)
    print("\nPrevalencia cruda por sistema (por 100,000) —", MAX_YEAR)
    print(crude)
    
    tabla_debug_sys = debug_crude_prevalence_by_system(base_prev, sys_pop_age, MAX_YEAR)
    print("\nResumen por previsión (numerador / denominador / tasa):")
    print(tabla_debug_sys)

    # ASR por sistema
    asr = who_asr_by_system(base_prev, sys_pop_age, MAX_YEAR)
    print("\nPrevalencia WHO ASR por sistema (por 100,000) —", MAX_YEAR)
    print(asr)

    # Conteos y tasas de prevalencia
    prev_obj = prevalence_counts(base_prev, pop_long)
    prev = prevalence_rates(prev_obj, pop_long)

    # Prevalencia cruda por año
    prev_crude = prev["crude_prevalence_per_year"].round(2).reindex(range(MIN_YEAR, MAX_YEAR + 1))
    prev_crude.to_csv(outdir / "crude_prevalence_per_year.csv", encoding="utf-8-sig")
    save_line(prev_crude, "Crude prevalence per year (per 100,000) – 2010–2023", gdir / "crude_prevalence_per_year.png")
    print("\nPREVALENCIA CRUDA (por 100,000) 2010–2023 [vivos = sin fecha fallecimiento]:")
    print(prev_crude.to_string())

    # WHO ASR prevalencia
    prev_asr = prev.get("who_asr_per_year", pd.Series(dtype=float)).reindex(range(MIN_YEAR, MAX_YEAR + 1))
    if prev_asr is not None and not prev_asr.dropna().empty:
        prev_asr.round(2).to_csv(outdir / "prevalence_who_asr_per_year.csv", encoding="utf-8-sig")
        save_line(prev_asr, "WHO-standardized prevalence (80+) (per 100,000) – 2010–2023",
                  gdir / "prevalence_who_asr_per_year.png")
        print("\nPREVALENCIA WHO-ESTANDARIZADA (80+) (por 100,000) 2010–2023:")
        print(prev_asr.round(2).to_string())

    # Gráfico comparativo incidencia vs prevalencia
    if not asr_inc.empty and not prev_asr.empty:
        save_two_series_line(
            asr_inc, prev_asr.loc[INC_START_YEAR:],
            ("Incidence (WHO ASR 80+)", "Prevalence (WHO ASR 80+)"),
            "Incidence vs Prevalence — WHO-standardized (per 100,000)",
            gdir / "incidence_prevalence_who_asr_comparative.png",
            x_index=range(INC_START_YEAR, MAX_YEAR + 1)
        )

    # ========================================
    # 5. ESTADÍSTICAS ADICIONALES
    # ========================================
    print("\n" + "=" * 60)
    print("ESTADÍSTICAS ADICIONALES...")
    print("=" * 60)
    
    # Reclasificación de previsión
    base_prev["prevision2"] = base_prev["prevision"].apply(
        lambda x: "ISAPRE" if x not in ("FONASA", "FFAA") else x
    )
    
    # Edad media por sistema
    stats_sistema = (
        base_prev
        .dropna(subset=["edad"])
        .groupby("prevision2")["edad"]
        .agg(n="count", mean="mean", sd="std")
    )
    
    stats_total = (
        base_prev["edad"]
        .dropna()
        .agg(n="count", mean="mean", sd="std")
        .to_frame().T
    )
    stats_total.index = ["TOTAL"]
    
    edad_stats_global = pd.concat([stats_sistema, stats_total]).round(2)
    print("\nEdad media y SD por sistema y total:")
    print(edad_stats_global)
    
    # Edad media por sexo
    stats_sexo = (
        base_prev
        .dropna(subset=["edad", "sexo"])
        .groupby("sexo")["edad"]
        .agg(n="count", mean="mean", sd="std")
    )
    
    stats_total_sexo = (
        base_prev["edad"]
        .dropna()
        .agg(n="count", mean="mean", sd="std")
        .to_frame().T
    )
    stats_total_sexo.index = ["TOTAL"]
    
    edad_stats_por_sexo = pd.concat([stats_sexo, stats_total_sexo]).round(2)
    print("\nEdad media y SD por sexo y total:")
    print(edad_stats_por_sexo)

    # ========================================
    # 6. GENERACIÓN DE GRÁFICOS
    # ========================================
    print("\n" + "=" * 60)
    print("GENERANDO GRÁFICOS...")
    print("=" * 60)
    
    # Pearson prevalencia vs latitud
    pearson_prevalence_vs_latitude(prev_obj, pop_long, base_prev, MAX_YEAR, outdir, gdir)

    # Histograma prevalencia por edad
    prevalence_age_hist_with_totals(base_prev, pop_long, MAX_YEAR, 
                                    gdir / "prevalence_age_hist_with_totals_2023.png")

    # Series temporales con IC
    plot_incidence_by_year_crude_and_asr_with_ci(incident, pop_long, gdir)
    plot_prevalence_by_year_crude_and_asr_with_ci(base_prev, pop_long, gdir)

    # Casos nuevos por sistema
    plot_new_cases_by_year_by_system(incident, gdir)

    # Detalle ISAPRE
    plot_isapre_detail_by_year(base_prev, gdir)

    # Edad promedio al diagnóstico
    plot_mean_age_at_incidence(incident, gdir)

    # Ratio y distribución por sexo
    plot_sex_ratio_and_percent_prevalent(prev_obj, gdir)

    # Prevalencia por sistema con IC
    plot_prevalence_by_system_over_time(base_prev, pop_long, gdir)

    # ========================================
    # 7. MAPAS REGIONALES
    # ========================================
    print("\n" + "=" * 60)
    print("GENERANDO MAPAS...")
    print("=" * 60)
    
    try:
        # Mapa TOTAL
        plot_prevalence_map_by_region(
            base_prev, pop_long, SHAPEFILE_PATH, MAX_YEAR,
            gdir / "map_prevalence_crude_by_region_TOTAL_2023.png",
            title="MS Prevalence by Region",
            filter_prevision=None,
            cmap_name="Reds"
        )

        # Mapa ISAPRE
        plot_prevalence_map_by_region(
            base_prev, pop_long, SHAPEFILE_PATH, MAX_YEAR,
            gdir / "map_prevalence_crude_by_region_ISAPRE_2023.png",
            title="MS Prevalence by Region",
            filter_prevision="ISAPRE",
            cmap_name="Greens",
            sys_pop_df=sys_pop_df
        )

        # Mapa FONASA
        plot_prevalence_map_by_region(
            base_prev, pop_long, SHAPEFILE_PATH, MAX_YEAR,
            gdir / "map_prevalence_crude_by_region_FONASA_2023.png",
            title="MS Prevalence by Region",
            filter_prevision="FONASA",
            cmap_name="Blues",
            sys_pop_df=sys_pop_df
        )
        
    except Exception as e:
        print(f"[WARN] No se pudieron generar mapas regionales: {e}")

    print("\n" + "=" * 60)
    print("ANÁLISIS COMPLETADO")
    print(f"Resultados guardados en: {outdir}")
    print(f"Gráficos guardados en: {gdir}")
    print("=" * 60)


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    main()

