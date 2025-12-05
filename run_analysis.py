#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Epidemiology of Relapsing-Remitting Multiple Sclerosis in Chile (2010-2023)
===========================================================================

Entry point script for running the complete epidemiological analysis.

Authors:
    - José M. Valdés (Clínica Alemana / UDD)
    - Javier Ortiz (FUVIEM)
    - Tamara Santibáñez (Clínica Alemana / UDD / Hospital La Florida)
    - Lorna Galleguillos (Clínica Alemana / UDD / Clínica Dávila)

Usage:
    python run_analysis.py

Requirements:
    - Python 3.10+
    - Dependencies: pandas, numpy, matplotlib, geopandas, scipy, statsmodels, unidecode, xlsxwriter, openpyxl

Input files required (in data/ subdirectories):
    - data/raw/BBDD EM.xlsx                                    # MS case database
    - data/population/estimaciones-y-proyecciones-*.xlsx       # INE projections
    - data/population/beneficiariosTramoEtario.xlsx            # Beneficiaries by age
    - data/population/Previsionporregion.xlsx                  # Population by insurance/region
    - data/geographic/Regional.shp (+ associated files)        # Regional shapefile

Outputs:
    - output/tables/    # CSV results
    - output/figures/   # PNG plots and maps
"""

import sys
from pathlib import Path

# Add src directory to path for imports
src_path = Path(__file__).parent / "src"
sys.path.insert(0, str(src_path))

from analisis_em.main import main

if __name__ == "__main__":
    main()
