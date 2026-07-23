# Dataset: data_MBA_modelling

## Overview
The **data_MBA_modelling** dataset contains 2,000 observations (rows) and 21 variables (columns). It is designed for statistical modelling or machine learning analyses, combining geographic, demographic, socio-economic, and serological data (antibody responses measured by Multiplex Bead Assay – MBA).

## Dataset Structure

| Feature          | Detail         |
|------------------|----------------|
| Number of rows   | 2,000          |
| Number of columns| 21             |
| Format           | Tabular (CSV/TSV) |

---

## Variable Dictionary

### Identifiers and Geographic/Environmental Data

| Variable        | Type      | Description |
|-----------------|-----------|-------------|
| `Units`         | Numeric   | Unique identifier for the statistical unit/individual. |
| `x`             | Numeric   | Spatial X coordinate (e.g., longitude or UTM coordinate). |
| `y`             | Numeric   | Spatial Y coordinate (e.g., latitude or UTM coordinate). |
| `pop.density`   | Numeric   | Population density in the area of residence. |
| `veg.stratum`   | Categorical | Vegetation stratum: `"Dense-vegetation"`, `"Moderate-vegetation"`, `"Sparse-vegetation"`. |

### Demographic and Socio-Economic Data

| Variable      | Type        | Description |
|---------------|-------------|-------------|
| `age`         | Numeric     | Age of the individual (in years). |
| `gender`      | Categorical | Gender: `"Male"`, `"Female"`. |
| `wealth`      | Categorical | Socio-economic status / wealth quintile: `"Wealthy"`, `"Middle-income"`, `"Low-income"`, `"Lower-middle-income"`. |
| `occupation`  | Categorical | Primary occupation: `"none"`, `"outside activities"`, `"housewife"`, `"other"`, `"student"`. |

### Serological Data – Multiplex Bead Assay (MBA)

The following variables represent antibody responses expressed as Median Fluorescence Intensity (MFI) against specific antigens from various pathogens.

| Variable                 | Pathogen / Antigen                                | Description |
|--------------------------|---------------------------------------------------|-------------|
| `Toxo.SAG2A`             | *Toxoplasma gondii* (SAG2A)                      | Antibody response to SAG2A antigen. |
| `Giardia.VSP3`           | *Giardia* (VSP3)                                 | Antibody response to VSP3 antigen. |
| `Giardia.VSP5`           | *Giardia* (VSP5)                                 | Antibody response to VSP5 antigen. |
| `Lf.Brugia.BmR1`         | Lymphatic filariasis – *Brugia malayi* (BmR1)   | Antibody response to BmR1 antigen. |
| `Lf.Wuchereria.Wb123`    | Lymphatic filariasis – *Wuchereria bancrofti* (Wb123) | Antibody response to Wb123 antigen. |
| `Lf.Brugia.Bm14`         | Lymphatic filariasis – *Brugia* (Bm14)           | Antibody response to Bm14 antigen. |
| `Lf.Brugia.Bm33`         | Lymphatic filariasis – *Brugia* (Bm33)           | Antibody response to Bm33 antigen. |
| `Trachoma.pgp3`          | Trachoma / *Chlamydia trachomatis* (pgp3)        | Antibody response to pgp3 antigen. |
| `Trachoma.ct694`         | Trachoma / *Chlamydia trachomatis* (ct694)       | Antibody response to ct694 antigen. |
| `Strongyloides.NIE`      | *Strongyloides stercoralis* (NIE)                | Antibody response to NIE antigen. |
| `Yaws.rp17`              | *Treponema pallidum pertenue* / Yaws (rp17)      | Antibody response to rp17 antigen. |
| `Yaws.TmpA`              | *Treponema pallidum pertenue* / Yaws (TmpA)      | Antibody response to TmpA antigen. |

---

## Usage Notes
- Serological variables are measured in **MFI** (Median Fluorescence Intensity) and can be used as predictors or outcomes in classification/regression models.
- Spatial coordinates (`x`, `y`) enable geospatial analysis or clustering.
- The `veg.stratum` variable is categorical, so it may require encoding for models that do not accept strings.
- This dataset is simulated for demonstration purposes; it does not represent real patient data.

---

## Example Loading in Python (pandas)
```python
import pandas as pd

df = pd.read_csv("data_MBA_modelling.csv")  # or .tsv
print(df.head())