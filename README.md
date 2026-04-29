# Barthel Index · Bayesian Clinical Decision Support

Interactive Shiny application for recovery probability estimation using a Bayesian hierarchical ordinal model on the Barthel Index of Activities of Daily Living.

---

## Branches

| Branch | Purpose |
|--------|---------|
| `main` | Stable, working code; production-ready |
| `experiments` | Active development: full MCMC settings (2000 iter, 4 chains), all features |
| `fast-experiments` | Lightweight version of `experiments`: fewer iterations (500), fewer chains (2), reduced `ndraws`; for quick feedback during development |

Changes flow: `fast-experiments` → `experiments` → `main` once validated.

---

## Pipelines

### Pipeline 1: Predictive Model (`01_model_training.R` + `app.R`)

Bayesian hierarchical ordinal regression predicting Barthel item recovery at discharge, conditioned on pathology, age, and comorbidity.

### Pipeline 2: Causal Inference (`02_causal_model_training.R`)

Estimates the **causal effect** of Cardiovascular Disease on Barthel outcomes using a Structural Causal Model (DAG), G-computation, and negative control outcome validation.

---

## Quick Start

### Prerequisites

- R >= 4.0
- Packages: `shiny`, `bslib`, `dplyr`, `ggplot2`, `brms`, `tidybayes`, `dagitty`, `marginaleffects`

```r
install.packages(c("shiny", "bslib", "dplyr", "ggplot2", "brms",
                   "tidybayes", "dagitty", "ggdag", "marginaleffects"))
```

### Step 1: Train the model

```bash
Rscript R/01_model_training.R
```

Runtime: Heavily depends on your machine (for a quick test try the version on branch `fastexperiments`). Outputs: `barthel_model.rds`, `session_info.txt`.

### Step 2: Launch the dashboard

```bash
Rscript -e "shiny::runApp('R')"
```

### Optional: Causal pipeline

```bash
Rscript R/02_causal_model_training.R
```

Outputs: `R/causal_barthel_model.rds`, `R/causal_dag.rds`, `R/causal_data_long.rds`.

---

## File Structure

```
R/
├── 01_model_training.R       # Predictive Bayesian model training
├── 02_causal_model_training.R # Causal inference pipeline (SCM + G-computation)
├── app.R                     # Shiny dashboard
data/
└── Dataset.csv               # Seccia et al. 2020 rehabilitation cohort
```

---

## Model Specification

```
Family:  cumulative(probit)
Formula: Score ~ Age + Pathology + Comorbidity + (1|PatientID) + (1|Item)
Items:   Feeding (max 10), Ambulation (max 15), Stairs (max 10)
```

Barthel items discretized into ordered categories: `Dependent < Assistance < Independent`.

---

## Troubleshooting

**`barthel_model.rds not found`**: Run `01_model_training.R` first.

**`Cannot locate data/Dataset.csv`**: Always run from the project root:
```bash
Rscript R/01_model_training.R   # correct
```

**Rhat > 1.01** — Increase `iter` and `chains` in the training script.
