#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Análisis Epidemiológico de Esclerosis Múltiple (EM) - Chile
============================================================

Script de entrada para ejecutar el análisis completo.

Uso:
    python run_analysis.py

Requisitos:
    - Python 3.10+
    - Dependencias: pandas, numpy, matplotlib, geopandas, scipy, statsmodels, unidecode, xlsxwriter, openpyxl

Archivos de entrada requeridos (en el mismo directorio):
    - BBDD EM.xlsx                                    # Base de datos de casos
    - estimaciones-y-proyecciones-2002-2035-comunas.xlsx  # Proyecciones INE
    - beneficiariosTramoEtario.xlsx                   # Beneficiarios por tramo etario
    - Previsionporregion.xlsx                         # Población por previsión y región
    - Regional.shp (+ archivos asociados)             # Shapefile de regiones

Salidas:
    - salidas_bdem/           # CSVs con resultados numéricos
    - salidas_bdem/graficos/  # Gráficos PNG
"""

from analisis_em import main

if __name__ == "__main__":
    main()

