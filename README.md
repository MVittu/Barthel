# Barthel Index · Bayesian Clinical Decision Support

> **Interactive Shiny Application for Recovery Probability Estimation**
>
> A Bayesian hierarchical ordinal model approach to overcome the "metric fallacy" of Barthel total scores by modeling individual recovery items with credible interval uncertainty quantification.

---

## Overview

This module implements a **clinical decision support system** for the Barthel Index of Activities of Daily Living (ADL). Rather than treating the Barthel total score as a continuous metric, we:

1. **Model individual items** (Feeding, Ambulation, Stairs) separately using ordinal regression
2. **Quantify prognostic uncertainty** via Bayesian 95% credible intervals (CrI)
3. **Provide recovery probabilities** stratified by patient pathology, age, and comorbidity
4. **Generate clinical alerts** based on item-level predictions

**Reference:** [Seccia et al. 2020 Dataset](https://figshare.com/articles/dataset/Data_of_patients_entering_a_rehabilitation_program/11663277) — stroke and orthopedic rehabilitation cohort

---

## Quick Start

### Prerequisites

- **R** ≥ 4.0
- **Packages:** `shiny`, `bslib`, `dplyr`, `ggplot2`, `brms`, `tidybayes`

Install all dependencies:

```r
packages <- c("shiny", "bslib", "dplyr", "ggplot2", "brms", "tidybayes")
install.packages(packages)
```

### Workflow

#### Step 1: Train the Bayesian Model

```bash
# From project root
Rscript R/01_model_training.R
```

Or in R/RStudio:
```r
source("R/01_model_training.R")
```

This script:
- Loads the raw dataset (`data/Dataset.csv`)
- Classifies ICD-9 diagnoses into **Pathology** (Neurological/Orthopedic) and **Comorbidity** (Cardiovascular/Diabetes/None)
- Transforms Barthel item scores into ordinal categories:
  - **Dependent:** Unable to perform activity
  - **Assistance:** Requires human or device assistance
  - **Independent:** Performs activity without assistance
- Fits a Bayesian cumulative-probit ordinal regression model with random effects for patients and items
- Saves the trained model as `barthel_model.rds` and session info as `session_info.txt`

**Runtime:** widely depends on your machine.

#### Step 2: Launch the Interactive Dashboard

```bash
# From project root
Rscript -e "shiny::runApp('R')"
```

Or in R/RStudio:
```r
setwd("R")
shiny::runApp()
```

The dashboard will open at `http://localhost:3838/` by default.

---

## Application Features

### Patient Profile Input

- **Main Pathology:** Neurological or Orthopedic
- **Age:** 30–90 years (slider)
- **Comorbidity:** None, Cardiovascular, or Diabetes

### Prediction Visualization

Three columns display **item-level recovery probabilities**:

- **Feeding** (blue): probability of independent feeding
- **Ambulation** (green): probability of independent ambulation
- **Stairs** (purple): probability of independent stair negotiation

Each chart shows:
- **Bar:** posterior mean probability
- **Error bars:** 95% credible interval

### Clinical Synthesis Panel

Biostatistical synthesis provides:
- **Overall prognosis** classification (Favourable/Moderate/Guarded)
- **Item-level probabilities** with 95% CrI and precision badges
- **Clinical alerts** for high uncertainty or poor outcomes:
  - Ambulation width >30pp: recommend re-evaluation
  - Stairs <25% probability: consider targeted intervention
  - Feeding uncertainty >30pp: OT assessment recommended

---

## File Structure

```
R/
├── README.md                 # This file
├── 01_model_training.R       # Bayesian model training pipeline
├── app.R                     # Shiny interactive dashboard
└── .gitkeep
```

### Generated Artifacts (after training)

```
. (project root)
├── barthel_model.rds         # Fitted brms model object
└── session_info.txt          # R/package versions for reproducibility
```

---

## Technical Details

### Model Specification

```
Family:    cumulative(probit)
Formula:   Score ~ Age + Pathology + Comorbidity + (1|PatientID) + (1|Item)
```

- **Cumulative-probit link:** Respects ordinal structure (Dependent < Assistance < Independent)
- **Random intercepts:**
  - `(1|PatientID)`: Patient-level heterogeneity
  - `(1|Item)`: Item difficulty (e.g., stairs harder than feeding)
- **Fixed effects:** Population-level trends by age, pathology type, and comorbidity

### Bayesian Inference Settings

- **Chains:** 2 (adaptive parallel sampling)
- **Iterations:** 2000 per chain
- **Warmup:** 1000 (discarded)
- **Posterior draws for prediction:** 500
- **Seed:** 42 (reproducible)

### Convergence Diagnostics

The training script checks **Rhat < 1.01** for all parameters to confirm MCMC convergence.

---

## Usage Examples

### Example 1: Neurological Patient, Age 72, No Comorbidity

1. Select **Pathology:** Neurological
2. Set **Age:** 72
3. Select **Comorbidity:** None
4. Click **Update Predictions**

**Expected output:** Posterior probabilities for Feeding, Ambulation, and Stairs, with clinical interpretation.

### Example 2: Orthopedic Patient Post-Hip Replacement

1. Select **Pathology:** Orthopedic
2. Set **Age:** 68
3. Select **Comorbidity:** Diabetes
4. Click **Update Predictions**

Dashboard will flag uncertainty in ambulation/stairs and recommend targeted assessments.

---

## Reproducibility & Dependencies

### Session Information

After model training, `session_info.txt` captures:
- R version
- OS and architecture
- All loaded packages and versions

**Note:** MCMC sampling is deterministic (`seed = 42`), but floating-point results may vary slightly across platforms due to BLAS/compiler differences. Statistical inference (credible intervals, ranks) remains robust.

### Verification

To confirm the model was fitted correctly:

```r
# In R console
barthel_model <- readRDS("barthel_model.rds")
summary(barthel_model)
# All Rhat < 1.01 confirms convergence
```

---

## Data Processing Details

### ICD-9 Classification Helpers

- **Pathology classification:** Prefix-based rules for stroke (430–438), CNS/PNS (320–359), musculoskeletal (710–739), fractures (800–829)
- **Comorbidity classification:** Scans secondary diagnoses for cardiovascular (401–448) or diabetes (250–259)

### Barthel Item Discretization

| Item | Dependent | Assistance | Independent |
|------|-----------|-----------|-------------|
| **Feeding** | Score = 0 | 0 < Score < 8 | Score ≥ 8 |
| **Ambulation** | Score ≤ 3 | 3 < Score < 12 | Score ≥ 12 |
| **Stairs** | Score ≤ 2 | 2 < Score < 8 | Score ≥ 8 |

---

## Troubleshooting

### Issue: `barthel_model.rds not found`

**Solution:** Run `01_model_training.R` first to generate the model.

```bash
Rscript R/01_model_training.R
```

### Issue: Rhat > 1.01 warnings

**Cause:** MCMC chain did not converge fully.  
**Solution:** Increase `iter` and `chains` in `01_model_training.R` (lines ~163–166).

### Issue: "Cannot locate data/Dataset.csv"

**Cause:** Script called from wrong working directory.  
**Solution:** Always run from the project root:

```bash
# Correct
Rscript R/01_model_training.R

# Wrong
cd R && Rscript 01_model_training.R
```

---

## References & Learning Resources

- **Bayesian Ordinal Regression:** [brms documentation](https://paul-buerkner.github.io/brms/)
- **Credible Intervals & Uncertainty:** [tidybayes cookbook](https://mjskay.github.io/tidybayes/)
- **Barthel Index:** WHO/FIM standardized ADL assessment
- **Seccia et al. 2020:** Reference dataset for stroke and orthopedic rehabilitation

---

## License & Attribution

This analysis builds on rehabilitation outcome research methodologies. Use responsibly for clinical audit and quality improvement.

**Citation:**  
Barthel Index · Bayesian Clinical Decision Support. [Work in progress].

---

## Questions?

For questions on model specification, Bayesian methods, or clinical interpretation, refer to:

- `01_model_training.R` comments (feature engineering, model setup)
- `app.R` comments (UI/server logic, posterior predictions)
- Inline code documentation throughout

---

**Last Updated:** 2026-04-28  
**Status:** Production Ready
