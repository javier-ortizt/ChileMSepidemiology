# Epidemiological Analysis of Multiple Sclerosis in Chile

Public code and aggregate outputs supporting the paper on relapsing-remitting multiple sclerosis (RRMS) in Chile, 2010–2023.

This public version is intentionally limited to:

- analysis code
- reproducible run scripts
- aggregate tables and figures used in the paper

It does **not** include raw or restricted source data.

## Scope

The repository supports the main epidemiological analyses reported in the paper:

- national incidence and prevalence
- WHO age-standardized rates
- public vs private insurance comparison
- latitude-gradient analysis
- age/sex-stratified analyses added during peer review

Key reported summary values in the current public outputs:

- Total unique cases: `4,118`
- National WHO-standardized prevalence in 2023: `19.77 per 100,000`
- Female proportion: `67.0%`
- ISAPRE prevalence in 2023: `61.60 per 100,000`
- FONASA prevalence in 2023: `12.97 per 100,000`
- Latitude correlation for WHO age-standardized prevalence: `r = 0.784`, `95% CI 0.471–0.921`, `p < 0.001`, `n = 16 regions`

## Public Repository Layout

```text
paper_epims/
├── README.md
├── LICENSE
├── CITATION.cff
├── requirements.txt
├── .gitignore
├── run_analysis.py
├── run_stratified_analysis.py
├── src/
│   └── analisis_em/
├── data/
│   └── README.md
├── results/
│   ├── tables/
│   ├── figures/
│   └── paper_assets/
└── output/                       # local generated artifacts, not tracked
```

## Restricted Inputs

The original case-level and denominator datasets are **not distributed** in this repository.

See [`data/README.md`](data/README.md) for the expected local file layout if you have lawful access to the underlying restricted inputs.

This means:

- the public repo is transparent about methods and code
- final aggregate outputs are included
- full reruns require local access to the non-public input files

## Code

Core analysis code lives in [`src/analisis_em`](src/analisis_em).

Main entry points:

- [`run_analysis.py`](run_analysis.py): full paper analysis pipeline
- [`run_stratified_analysis.py`](run_stratified_analysis.py): stratified age/sex analysis used for reviewer-requested additions

## Public Results

Tracked paper-facing outputs are stored in [`results/`](results).

### Tables

Files in [`results/tables`](results/tables):

- `crude_incidence_per_year.csv`
- `incidence_who_asr_per_year.csv`
- `crude_prevalence_per_year.csv`
- `prevalence_who_asr_per_year.csv`
- `prevalence_latitude_merge_2023.csv`
- `table_incidence_by_age_sex.csv`
- `table_prevalence_by_age_sex_2023.csv`
- `table_regional_asr_prevalence_2023.csv`

### Figures

Files in [`results/figures`](results/figures):

- `crude_incidence_per_year.png`
- `crude_prevalence_per_year.png`
- `map_prevalence_crude_by_region_TOTAL_2023.png`
- `prevalence_pearson_who_asr_vs_latitude_2023.png`
- `fig_incidence_by_age_sex.png`
- `fig_prevalence_by_age_sex_2023.png`

### Paper Assets

Files in [`results/paper_assets`](results/paper_assets):

- `Tables_4_5_6.docx`

## Mapping To Paper Elements

The tracked public outputs correspond to the paper as follows:

- Main incidence time series: `results/tables/crude_incidence_per_year.csv` and `results/figures/crude_incidence_per_year.png`
- Main prevalence time series: `results/tables/crude_prevalence_per_year.csv` and `results/figures/crude_prevalence_per_year.png`
- Latitude/regional analysis: `results/tables/prevalence_latitude_merge_2023.csv`, `results/tables/table_regional_asr_prevalence_2023.csv`, `results/figures/map_prevalence_crude_by_region_TOTAL_2023.png`, and `results/figures/prevalence_pearson_who_asr_vs_latitude_2023.png`
- Reviewer-added age/sex stratification: `results/tables/table_incidence_by_age_sex.csv`, `results/tables/table_prevalence_by_age_sex_2023.csv`, `results/figures/fig_incidence_by_age_sex.png`, and `results/figures/fig_prevalence_by_age_sex_2023.png`

## Local Reproduction

Install dependencies:

```bash
py -m venv .venv
.venv\Scripts\activate
py -m pip install -r requirements.txt
```

Run the main analysis:

```bash
py run_analysis.py
```

Run the stratified reviewer-response analysis:

```bash
py run_stratified_analysis.py
```

Both scripts expect the restricted local inputs described in [`data/README.md`](data/README.md).

## Methods Summary

- Incident case: first eligible record per beneficiary
- Prevalent case: diagnosed up to year `Y` and alive at end of `Y`
- Incidence baseline handling: FFAA cases before 2010 are consolidated into 2010; incidence reporting excludes 2010 as baseline where applicable
- Age standardization: WHO 2000–2025 World Standard Population with open `80+` group
- Latitude analysis: Pearson correlation with Fisher z confidence interval

## License

This project is licensed under the MIT License. See [`LICENSE`](LICENSE).

## Funding

This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.
