# Epidemiology of Relapsing-Remitting Multiple Sclerosis in Chile (2010–2023)

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the analysis code, data, and methodology for the epidemiological study of **Relapsing-Remitting Multiple Sclerosis (RRMS)** in Chile using administrative health records from the GES (Explicit Health Guarantees) system.

## Authors

- **José M. Valdés**<sup>1,2</sup> ✉️ (josemiguelvv@gmail.com)
- **Javier Ortiz**<sup>3</sup>
- **Tamara Santibáñez**<sup>1,2,4</sup>
- **Lorna Galleguillos**<sup>1,2,5</sup>

### Affiliations

1. Clínica Alemana de Santiago, Department of Neurology and Psychiatry, Santiago, Chile
2. Facultad de Medicina, Clínica Alemana Universidad del Desarrollo, Santiago, Chile
3. Fundación para la Vida con Esclerosis Múltiple (FUVIEM), Chile
4. Hospital Clínico Dra. Eloísa Díaz Insunza, La Florida, Chile
5. Clínica Dávila, Neurology and Neurosurgery Department, Santiago, Chile

---

## Table of Contents

- [Abstract](#abstract)
- [Key Findings](#key-findings)
- [Repository Structure](#repository-structure)
- [Data Sources](#data-sources)
- [Methodology](#methodology)
- [Installation](#installation)
- [Usage](#usage)
- [Output Files](#output-files)
- [Statistical Methods](#statistical-methods)
- [Limitations](#limitations)
- [Citation](#citation)
- [License](#license)

---

## Abstract

**Background:** Multiple sclerosis (MS) is a chronic neuroinflammatory disease, yet comprehensive, nationwide epidemiological data from Chile remain scarce.

**Objectives:** To estimate the national incidence and prevalence of MS in Chile from 2010 to 2023, characterize demographic patterns, and quantify disparities across the nation's segmented public and private healthcare systems.

**Methods:** A nationwide, retrospective analysis of administrative health data was conducted. MS cases diagnosed from 2010 to 2023 were identified from the databases of the National Health Fund (FONASA), private insurance institutions (ISAPREs), and the Armed Forces Health System, covering the entire population. Incidence and prevalence rates were calculated using population projections from the National Statistical Institute (INE) and were age-standardized to the WHO 2000–2025 World Standard Population.

**Results:** A total of 4,118 unique MS cases were identified. The national age-adjusted prevalence of MS in 2023 was **19.77 per 100,000** inhabitants. The cohort was predominantly female (67.0%), with a female-to-male ratio of 2.03:1. The mean age at diagnosis was 35.99 years (SD ±12.6). A striking disparity in diagnosed prevalence was observed between the private ISAPRE system (61.60 per 100,000) and the public FONASA system (12.97 per 100,000), representing a **4.75-fold difference**. Furthermore, a significant positive linear relationship between MS prevalence and increasing southern latitude was confirmed (Pearson's r = 0.766, p < 0.001).

**Conclusion:** This study provides the first comprehensive epidemiological profile of MS in Chile, revealing a moderate national prevalence and confirming a latitudinal gradient. The most critical finding is a profound disparity in diagnosed prevalence between the private and public healthcare systems.

**Keywords:** Multiple Sclerosis, Chile, Epidemiology, Healthcare Disparities, Latitude Gradient

---

## Key Findings

| Metric | Value |
|--------|-------|
| Total unique cases | 4,118 |
| National WHO-adjusted prevalence (2023) | 19.77 per 100,000 |
| Female proportion | 67.0% |
| Female-to-male ratio | 2.03:1 |
| Mean age at diagnosis | 35.99 years (SD ±12.6) |
| ISAPRE prevalence (2023) | 61.60 per 100,000 |
| FONASA prevalence (2023) | 12.97 per 100,000 |
| ISAPRE/FONASA ratio | 4.75× |
| Latitude correlation (Pearson r) | 0.766 (p < 0.001) |

---

## Repository Structure

```
ChileMSepidemiology/
│
├── README.md                       # This file
├── LICENSE                         # MIT License
├── requirements.txt                # Python dependencies
├── .gitignore                      # Git ignore rules
│
├── src/                            # Source code
│   └── analisis_em/                # Main analysis package
│       ├── __init__.py             # Package initialization
│       ├── config.py               # Configuration and constants
│       ├── utils.py                # Utility functions
│       ├── io_data.py              # Data loading functions
│       ├── who_standards.py        # WHO standardization weights
│       ├── geography.py            # Geographic helpers
│       ├── statistics.py           # Statistical methods
│       ├── incidence.py            # Incidence calculations
│       ├── prevalence.py           # Prevalence calculations
│       ├── insurance.py            # Insurance system analysis
│       ├── correlation.py          # Latitude correlation
│       ├── plotting.py             # Visualization functions
│       ├── maps.py                 # Choropleth map generation
│       └── main.py                 # Main orchestration
│
├── data/                           # Input data files
│   ├── raw/                        # Original unprocessed data
│   │   ├── BBDD EM.xlsx            # MS case registry (GES)
│   │   └── fonasa2010-2022.xlsx    # Historical FONASA data
│   │
│   ├── population/                 # Population denominators
│   │   ├── estimaciones-y-proyecciones-2002-2035-comunas.xlsx  # INE projections
│   │   ├── beneficiariosTramoEtario.xlsx   # Beneficiaries by age group
│   │   └── Previsionporregion.xlsx         # Population by insurance/region
│   │
│   └── geographic/                 # Geographic data
│       ├── Regional.shp            # Chile regions shapefile
│       ├── Regional.shx            # Shapefile index
│       ├── Regional.dbf            # Shapefile attributes
│       ├── Regional.prj            # Projection file
│       └── coordenadas.xlsx        # Regional coordinates
│
├── output/                         # Generated results (gitignored)
│   ├── tables/                     # CSV results
│   └── figures/                    # Generated plots and maps
│
├── docs/                           # Documentation
│   └── methodology.md              # Detailed methodology
│
└── run_analysis.py                 # Entry point script
```

---

## Data Sources

| Dataset | Description | Source |
|---------|-------------|--------|
| `BBDD EM.xlsx` | MS case registry from GES system (2010–2023) | FONASA, ISAPREs, FFAA |
| `estimaciones-y-proyecciones-2002-2035-comunas.xlsx` | Population projections by age, sex, and region | National Statistics Institute (INE) |
| `beneficiariosTramoEtario.xlsx` | Health insurance beneficiaries by age group | FONASA/ISAPRE registries |
| `Previsionporregion.xlsx` | Population by insurance system and region | Health authorities |
| `Regional.shp` | Administrative boundaries shapefile (16 regions) | Chilean geographic data |

### Case Definition

- **Incident case**: First recorded entry in the GES system per unique beneficiary
- **Prevalent case**: All cases diagnosed up to reference year who are **alive** (no recorded death date)
- **FFAA adjustment**: Cases from Armed Forces system prior to 2010 are reassigned to 2010 baseline

### Healthcare System Context

In 2010, Relapsing-Remitting MS was incorporated into Chile's GES (Garantías Explícitas en Salud) system, mandating the reporting of all diagnosed cases. Chile's healthcare system is segmented into:

- **FONASA** (public): ~80% of population
- **ISAPREs** (private): ~17% of population  
- **FFAA** (armed forces): ~3% of population

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
- **Denominator**: INE population projections (total) or insurance-specific population
- **Period**: 2010–2023

### Age Standardization

WHO World Standard Population weights (2000–2025) with 80+ collapsed group:

| Age Group | Weight |
|-----------|--------|
| 0–4       | 0.0886 |
| 5–9       | 0.0869 |
| 10–14     | 0.0860 |
| 15–19     | 0.0847 |
| 20–24     | 0.0822 |
| 25–29     | 0.0793 |
| 30–34     | 0.0761 |
| 35–39     | 0.0715 |
| 40–44     | 0.0659 |
| 45–49     | 0.0604 |
| 50–54     | 0.0537 |
| 55–59     | 0.0455 |
| 60–64     | 0.0372 |
| 65–69     | 0.0296 |
| 70–74     | 0.0221 |
| 75–79     | 0.0152 |
| 80+       | 0.0154 |

```
ASR = Σ (age-specific rate × WHO weight) / Σ (WHO weights for available groups)
```

---

## Installation

### Requirements

- Python 3.10 or higher
- Dependencies listed in `requirements.txt`

### Setup

```bash
# Clone the repository
git clone https://github.com/javier-ortizt/ChileMSepidemiology.git
cd ChileMSepidemiology

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
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
4. WHO age standardization
5. Insurance system stratification (FONASA vs ISAPRE)
6. Geographic/latitude analysis
7. Plot and map generation

### Using Individual Modules

```python
from src.analisis_em.io_data import read_bdem, load_ine_population
from src.analisis_em.incidence import compute_incidence_from_bdem
from src.analisis_em.prevalence import build_prevalence_base, prevalence_counts

# Load data
df = read_bdem("data/raw/BBDD EM.xlsx")
pop_long = load_ine_population("data/population/estimaciones-y-proyecciones-2002-2035-comunas.xlsx")

# Calculate incidence
incident = compute_incidence_from_bdem(df)

# Calculate prevalence
base_prev = build_prevalence_base(df)
prev_obj = prevalence_counts(base_prev, pop_long)
```

---

## Output Files

### Tables (`output/tables/`)

| File | Description |
|------|-------------|
| `crude_incidence_per_year.csv` | Crude incidence rates per 100,000 (2011–2023) |
| `incidence_who_asr_per_year.csv` | WHO age-standardized incidence |
| `crude_prevalence_per_year.csv` | Crude prevalence rates per 100,000 (2010–2023) |
| `prevalence_who_asr_per_year.csv` | WHO age-standardized prevalence |
| `prevalence_by_system.csv` | Prevalence by FONASA/ISAPRE with 95% CI |
| `prevalence_latitude_correlation.csv` | Regional prevalence with latitude data |

### Figures (`output/figures/`)

| File | Description |
|------|-------------|
| `incidence_crude_by_sex_total_with_ci.png` | Crude incidence time series with 95% CI |
| `prevalence_crude_by_sex_total_with_ci.png` | Crude prevalence time series with 95% CI |
| `prevalence_who_asr_by_sex_total_with_ci.png` | ASR prevalence with 95% CI |
| `prevalence_pearson_crude_vs_latitude_2023.png` | Latitude correlation scatter plot |
| `map_prevalence_crude_by_region_*.png` | Choropleth maps by region |
| `prevalence_by_system_comparison.png` | FONASA vs ISAPRE comparison |

---

## Statistical Methods

### Confidence Intervals

**Crude rates (Poisson):**
```
Lower = 0.5 × χ²(α/2, 2×count)
Upper = 0.5 × χ²(1-α/2, 2×(count+1))
```

**ASR variance (approximate):**
```
Var(ASR) = Σᵢ wᵢ² × (dᵢ / Nᵢ²) × 100,000²
95% CI = ASR ± 1.96 × √Var(ASR)
```

**Pearson correlation (Fisher's z):**
```
z = arctanh(r)
SE(z) = 1/√(n-3)
CI(r) = tanh(z ± 1.96 × SE(z))
```

---

## Limitations

1. **Administrative data**: Cases identified through GES system; may underestimate true prevalence in populations with limited healthcare access
2. **Survival assumption**: Cases without death date assumed alive; may slightly overestimate prevalence
3. **FFAA historical data**: Pre-2010 Armed Forces cases reassigned to 2010 baseline
4. **Age at diagnosis**: Age recorded at GES entry, not necessarily at symptom onset
5. **Diagnostic disparities**: The 4.75× difference between ISAPRE and FONASA likely reflects access inequities, not true disease prevalence differences

---

## Reproducibility

All analysis parameters are centralized in `src/analisis_em/config.py`:

```python
MIN_YEAR = 2010          # Start year for prevalence
INC_START_YEAR = 2010    # Baseline year for incidence (excluded from reporting)
MAX_YEAR = 2023          # End year
TOP_OPEN = 80            # Open age group (80+)
```

---

## Citation

If you use this code or data, please cite:

```bibtex
@article{valdes2025ms_chile,
  title={Epidemiology of Relapsing Remitting Multiple Sclerosis in Chile},
  author={Valdés, José M. and Ortiz, Javier and Santibáñez, Tamara and Galleguillos, Lorna},
  journal={[Journal]},
  year={2025},
  note={In preparation}
}
```

---

## Funding

This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contact

For questions regarding the analysis methodology or code:

**Corresponding Author:**  
Dr. José M. Valdés  
📧 josemiguelvv@gmail.com  
Neurology and Psychiatry Department, Clínica Alemana  
Av. Vitacura 5951, 10th Floor, Vitacura, Santiago, Chile

---

*Last updated: December 202*
