# Epidemiological Analysis of Multiple Sclerosis in Chile (2010–2023)

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the analysis code and methodology for calculating **incidence and prevalence rates of Multiple Sclerosis (MS)** in Chile using administrative health records from the GES (Explicit Health Guarantees) system.

## Table of Contents

- [Overview](#overview)
- [Data Sources](#data-sources)
- [Methodology](#methodology)
- [Installation](#installation)
- [Usage](#usage)
- [Output Files](#output-files)
- [Project Structure](#project-structure)
- [Statistical Methods](#statistical-methods)
- [Limitations](#limitations)
- [Citation](#citation)

---

## Overview

This study estimates the epidemiological burden of Multiple Sclerosis in Chile by calculating:

- **Crude incidence rates** (2011–2023)
- **Crude prevalence rates** (2010–2023)
- **Age-standardized rates (ASR)** using the WHO World Standard Population (2000–2025)
- **Regional prevalence** with latitude correlation analysis
- **Prevalence by health insurance system** (FONASA vs ISAPRE)

### Key Features

- Age-specific rates by 5-year age groups (0–4, 5–9, ..., 75–79, 80+)
- Sex-stratified analysis (Male/Female)
- 95% confidence intervals using exact Poisson methods
- WHO age-standardization with collapsed 80+ group
- Geographic analysis by Chile's 16 administrative regions
- Pearson correlation between prevalence and geographic latitude

---

## Data Sources

| Dataset | Description | Source |
|---------|-------------|--------|
| `BBDD EM.xlsx` | MS case registry from GES system | Administrative health records |
| `estimaciones-y-proyecciones-2002-2035-comunas.xlsx` | Population projections by age, sex, and region | National Statistics Institute (INE) |
| `beneficiariosTramoEtario.xlsx` | Health insurance beneficiaries by age group | FONASA/ISAPRE registries |
| `Previsionporregion.xlsx` | Population by insurance system and region | Health authorities |
| `Regional.shp` | Administrative boundaries shapefile | Geographic data portal |

### Case Definition

- **Incident case**: First recorded entry in the GES system per unique beneficiary
- **Prevalent case**: All cases diagnosed up to reference year who are **alive** (no recorded death date)
- **FFAA adjustment**: Cases from Armed Forces (FFAA) system prior to 2010 are reassigned to 2010 baseline

---

## Methodology

### Incidence Calculation

```
Crude Incidence Rate = (New cases in year Y / Mid-year population) × 100,000
```

- **Numerator**: First GES registration per beneficiary
- **Denominator**: INE population projections
- **Period**: 2011–2023 (2010 used as baseline, excluded from reporting)

### Prevalence Calculation

```
Crude Prevalence Rate = (Alive cases diagnosed ≤ year Y / Mid-year population) × 100,000
```

- **Numerator**: Cumulative cases without death date recorded
- **Denominator**: INE population projections
- **Period**: 2010–2023

### Age Standardization

WHO World Standard Population weights (2000–2025) with 80+ collapsed group:

| Age Group | Weight |
|-----------|--------|
| 0–4       | 0.0886 |
| 5–9       | 0.0869 |
| ...       | ...    |
| 75–79     | 0.0152 |
| 80+       | 0.0154 |

```
ASR = Σ (age-specific rate × WHO weight) / Σ (WHO weights for available groups)
```

### Confidence Intervals

- **Crude rates**: Exact Poisson CI using chi-square limits
- **ASR**: Approximate variance method: `Var(ASR) ≈ Σ wᵢ² × (dᵢ / Nᵢ²) × 100,000²`
- **Pearson correlation**: Fisher's z-transformation

---

## Installation

### Requirements

- Python 3.10 or higher
- Required packages listed in `requirements.txt`

### Setup

```bash
# Clone or download the repository
cd paper_epims

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Dependencies

```
pandas>=2.0.0
numpy>=1.24.0
matplotlib>=3.7.0
geopandas>=0.13.0
scipy>=1.10.0
statsmodels>=0.14.0
openpyxl>=3.1.0
xlsxwriter>=3.1.0
unidecode>=1.3.0
```

---

## Usage

### Run Complete Analysis

```bash
python run_analysis.py
```

This executes the full analysis pipeline:
1. Data loading and validation
2. Incidence calculation (2011–2023)
3. Prevalence calculation (2010–2023)
4. Age standardization (WHO ASR)
5. Insurance system stratification
6. Geographic analysis
7. Plot and map generation

### Using Individual Modules

```python
from analisis_em.io_data import read_bdem, load_ine_population
from analisis_em.incidence import compute_incidence_from_bdem, incidence_rates_and_counts
from analisis_em.prevalence import build_prevalence_base, prevalence_counts, prevalence_rates
from analisis_em.statistics import poisson_rate_ci, asr_and_ci_from_counts

# Load data
df = read_bdem("BBDD EM.xlsx")
pop_long = load_ine_population("estimaciones-y-proyecciones-2002-2035-comunas.xlsx")

# Calculate incidence
incident = compute_incidence_from_bdem(df)
inc_rates = incidence_rates_and_counts(incident, pop_long)

# Calculate prevalence
base_prev = build_prevalence_base(df)
prev_obj = prevalence_counts(base_prev, pop_long)
prev_rates = prevalence_rates(prev_obj, pop_long)
```

---

## Output Files

### CSV Results (`salidas_bdem/`)

| File | Description |
|------|-------------|
| `new_cases_per_year.csv` | Annual incident case counts |
| `crude_incidence_per_year.csv` | Crude incidence rates per 100,000 |
| `incidence_who_asr_per_year.csv` | WHO age-standardized incidence |
| `age_specific_incidence.csv` | Incidence by year and age group |
| `crude_prevalence_per_year.csv` | Crude prevalence rates per 100,000 |
| `prevalence_who_asr_per_year.csv` | WHO age-standardized prevalence |
| `prevalence_latitude_merge_2023.csv` | Regional prevalence with latitude |

### Figures (`salidas_bdem/graficos/`)

| File | Description |
|------|-------------|
| `incidence_crude_by_sex_total_with_ci.png` | Crude incidence time series with 95% CI |
| `incidence_who_asr_by_sex_total_with_ci.png` | ASR incidence with 95% CI |
| `prevalence_crude_by_sex_total_with_ci.png` | Crude prevalence time series with 95% CI |
| `prevalence_who_asr_by_sex_total_with_ci.png` | ASR prevalence with 95% CI |
| `prevalence_pearson_crude_vs_latitude_2023.png` | Latitude correlation (crude) |
| `prevalence_pearson_who_asr_vs_latitude_2023.png` | Latitude correlation (ASR) |
| `map_prevalence_crude_by_region_*.png` | Choropleth maps by region |
| `new_cases_by_system_*.png` | Cases by insurance system |
| `mean_age_at_incidence_by_year.png` | Mean age at diagnosis |
| `sex_ratio_m_over_f_prevalent_by_year.png` | Sex ratio trends |

---

## Project Structure

```
paper_epims/
│
├── analisis_em/                    # Main analysis package
│   ├── __init__.py                 # Package initialization
│   ├── config.py                   # Configuration and constants
│   ├── utils.py                    # Utility functions
│   ├── io_data.py                  # Data loading functions
│   ├── who_standards.py            # WHO standardization weights
│   ├── geography.py                # Geographic helpers
│   ├── statistics.py               # Statistical methods
│   ├── incidence.py                # Incidence calculations
│   ├── prevalence.py               # Prevalence calculations
│   ├── insurance.py                # Insurance system analysis
│   ├── correlation.py              # Latitude correlation
│   ├── plotting.py                 # Visualization functions
│   ├── maps.py                     # Choropleth map generation
│   └── main.py                     # Main orchestration
│
├── run_analysis.py                 # Entry point script
├── requirements.txt                # Python dependencies
├── README.md                       # This file
│
├── BBDD EM.xlsx                    # Input: MS case database
├── estimaciones-y-proyecciones-2002-2035-comunas.xlsx  # Input: INE population
├── beneficiariosTramoEtario.xlsx   # Input: Beneficiaries by age
├── Previsionporregion.xlsx         # Input: Population by insurance
├── Regional.shp                    # Input: Regional shapefile
│
└── salidas_bdem/                   # Output directory
    ├── *.csv                       # Numeric results
    └── graficos/                   # Generated figures
        └── *.png
```

---

## Statistical Methods

### Poisson Confidence Intervals

For count data, exact Poisson confidence intervals are calculated using chi-square quantiles:

```
Lower = 0.5 × χ²(α/2, 2×count)
Upper = 0.5 × χ²(1-α/2, 2×(count+1))
```

### Age-Standardized Rate (ASR) Variance

```
Var(ASR) = Σᵢ wᵢ² × (dᵢ / Nᵢ²) × multiplier²
SE(ASR) = √Var(ASR)
95% CI = ASR ± 1.96 × SE(ASR)
```

Where:
- wᵢ = WHO weight for age group i (re-normalized)
- dᵢ = cases in age group i
- Nᵢ = population in age group i

### Pearson Correlation CI

Fisher's z-transformation:

```
z = arctanh(r)
SE(z) = 1/√(n-3)
CI(z) = z ± z_α/2 × SE(z)
CI(r) = tanh(CI(z))
```

---

## Limitations

1. **Administrative data**: Cases identified through GES system; may underestimate true prevalence
2. **Survival assumption**: Cases without death date assumed alive; may overestimate prevalence
3. **FFAA historical data**: Pre-2010 Armed Forces cases reassigned to 2010 baseline
4. **Age at diagnosis**: Age recorded at GES entry, not necessarily at symptom onset
5. **Regional population**: Insurance-specific denominators limited to available data

---

## Reproducibility

To ensure reproducibility:

1. All random seeds are fixed where applicable
2. Data processing steps are deterministic
3. Version requirements specified in `requirements.txt`
4. Analysis parameters centralized in `config.py`

### Key Parameters

```python
MIN_YEAR = 2010          # Start year for prevalence
INC_START_YEAR = 2010    # Baseline year for incidence (excluded from reporting)
MAX_YEAR = 2023          # End year
TOP_OPEN = 80            # Open age group (80+)
```

---


---

## Contact

For questions regarding the analysis methodology or code, please contact the corresponding author.

---

## License

This project is licensed under the MIT License - see the LICENSE file for details.

---

*Last updated: December 2024*

