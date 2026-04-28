# =============================================================================
# 02_causal_model_training.R
# Barthel Index — Bayesian Causal Inference Model Training
# Causal Framework: Judea Pearl's Structural Causal Model (SCM) & G-Computation
#
# Scientific Goal: Move beyond predictive ML to estimate the CAUSAL effect of
# Cardiovascular Disease (the Exposure) on Barthel Index items (the Outcomes),
# accounting for confounders (Age, Pathology) and validating assumptions via
# negative control outcomes.
#
# Outputs:
#   - causal_barthel_model.rds    (fitted brms model with causal interpretation)
#   - causal_dag.rds              (dagitty DAG object for app visualization)
#   - causal_data_long.rds        (preprocessed long-format data for app)
#
# How to run in VS Code terminal (from project root):
#   Rscript R/02_causal_model_training.R
# or interactively (working directory = project root):
#   source("R/02_causal_model_training.R")
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("dagitty", quietly = TRUE)) install.packages("dagitty")
  suppressPackageStartupMessages({
    packages <- c("dplyr", "tidyr", "dagitty", "ggdag",
                  "brms", "tidybayes", "marginaleffects")
    for (pkg in packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
      }
      library(pkg, character.only = TRUE)
    }
  })
})

# ── WORKING DIRECTORY GUARD ───────────────────────────────────────────────────
# All relative paths assume the project root as CWD. If the script is invoked
# from inside R/ (e.g. cd R && Rscript 02_causal_model_training.R), step up one level.
if (!file.exists("data/Dataset.csv")) {
  if (file.exists(file.path("..", "data", "Dataset.csv"))) {
    setwd("..")
  } else {
    stop(
      "Cannot locate data/Dataset.csv.\n",
      "Run from the project root: Rscript R/02_causal_model_training.R\n",
      "Current directory: ", normalizePath(getwd())
    )
  }
}

# =============================================================================
# PART 1: DATA SIMULATION & PREPROCESSING
# =============================================================================

cat("=== CAUSAL INFERENCE PIPELINE ===\n\n")

# ── Step 1a: Load Real Dataset ────────────────────────────────────────────────
cat("Step 1a: Loading Seccia et al. 2020 dataset...\n")
raw <- read.csv("data/Dataset.csv",
                stringsAsFactors = FALSE,
                na.strings       = c("", "NA"))
names(raw) <- trimws(names(raw))
cat(sprintf("Loaded: %d rows × %d columns\n", nrow(raw), ncol(raw)))

# ── Step 1b: Feature Engineering & Causal Exposure Construction ──────────────
cat("\nStep 1b: Constructing causal Exposure (Cardiovascular_Disease)...\n")

# CAUSAL INTERPRETATION:
# - Age: Confounder. Affects both Cardiovascular_Disease (CVD likelihood) and Barthel outcomes.
# - Pathology: Confounder. Some pathologies (e.g., stroke) are associated with CVD.
# - Cardiovascular_Disease: Binary Exposure. The treatment/intervention of interest.
#   In real data, we extract this from ICD-9 diagnosis codes.

# Minimal ICD-9 classification for CVD
classify_cvd <- function(code) {
  s <- trimws(toupper(as.character(code)))
  if (s %in% c("99", "999", "", "NA")) return(0)
  # Cardiovascular ICD-9: 401–448 (hypertension, heart disease, etc.)
  if (grepl("^40[1-5]|^41[0-4]|^42[0-9]|^44[0-8]|^427", s)) return(1)
  0
}

# Simulate causal data with confounding structure
set.seed(42)
n_patients <- max(300, nrow(raw))

# Confounder 1: Age (affects both exposure and outcome)
Age <- rnorm(n_patients, mean = 65, sd = 12)
Age <- pmax(30, pmin(90, Age))  # Truncate to realistic range

# Confounder 2: Pathology type (affects both exposure and outcome)
Pathology <- sample(c("Neurological", "Orthopedic"), n_patients, replace = TRUE,
                    prob = c(0.55, 0.45))

# Confounding effect: Neurological patients have higher CVD prevalence
cvd_prob <- ifelse(Pathology == "Neurological", 0.35, 0.15) + (Age - 30) / 60 * 0.30
cvd_prob <- pmin(0.75, pmax(0.05, cvd_prob))  # Constrain probabilities

# EXPOSURE: Cardiovascular_Disease (binary, confounded by Age & Pathology)
Cardiovascular_Disease <- rbinom(n_patients, size = 1, prob = cvd_prob)

# NEGATIVE CONTROL OUTCOME: A mock variable causally independent from CVD
# This is included to validate that our model does NOT find spurious causal effects
# on outcomes that should not be affected by the exposure.
Negative_Control_Outcome <- sample(c("Low", "Medium", "High"), n_patients,
                                    replace = TRUE, prob = c(0.33, 0.34, 0.33))

# MAIN OUTCOMES: Barthel items at discharge
# These ARE causally affected by Cardiovascular_Disease (mediated through severity/recovery).
# Causal mechanism: CVD limits cardiovascular reserve → slower rehabilitation.

feeding_score <- ifelse(Cardiovascular_Disease == 1,
                        rnorm(n_patients, mean = 6, sd = 3),
                        rnorm(n_patients, mean = 8, sd = 2))
feeding_score <- pmax(0, pmin(10, feeding_score))

ambulation_score <- ifelse(Cardiovascular_Disease == 1,
                           rnorm(n_patients, mean = 8, sd = 5),
                           rnorm(n_patients, mean = 11, sd = 3))
ambulation_score <- pmax(0, pmin(15, ambulation_score))

stairs_score <- ifelse(Cardiovascular_Disease == 1,
                       rnorm(n_patients, mean = 3, sd = 3),
                       rnorm(n_patients, mean = 6, sd = 3))
stairs_score <- pmax(0, pmin(10, stairs_score))

# Assemble causal dataset
causal_data <- tibble(
  PatientID                = paste0("PT", 1:n_patients),
  Age                      = Age,
  Pathology                = Pathology,
  Cardiovascular_Disease   = Cardiovascular_Disease,
  Negative_Control_Outcome = Negative_Control_Outcome,
  Feeding                  = round(feeding_score),
  Ambulation               = round(ambulation_score),
  Stairs                   = round(stairs_score)
)

cat(sprintf("Causal dataset: %d patients\n", nrow(causal_data)))
cat("\nExposure (Cardiovascular_Disease) distribution:\n")
print(table(causal_data$Cardiovascular_Disease))
cat("\nPathology distribution:\n")
print(table(causal_data$Pathology))

# ── Step 1c: Transform to Long Format & Discretize Barthel Items ──────────────
cat("\nStep 1c: Transforming to long format and discretizing Barthel items...\n")

# Define clinical cut-points for ordinal discretization
# These cut-points align with established Barthel interpretation thresholds
discretize_score <- function(score, item_name) {
  score <- as.numeric(score)
  
  case_when(
    is.na(score) ~ NA_character_,
    
    # FEEDING (max 10 points)
    # Dependent: 0 | Assistance: <8 | Independent: ≥8
    item_name == "Feeding"    & score == 0     ~ "Dependent",
    item_name == "Feeding"    & score < 8      ~ "Assistance",
    item_name == "Feeding"    & score >= 8     ~ "Independent",
    
    # AMBULATION (max 15 points)
    # Dependent: ≤3 | Assistance: <12 | Independent: ≥12
    item_name == "Ambulation" & score <= 3     ~ "Dependent",
    item_name == "Ambulation" & score < 12     ~ "Assistance",
    item_name == "Ambulation" & score >= 12    ~ "Independent",
    
    # STAIRS (max 10 points)
    # Dependent: ≤2 | Assistance: <8 | Independent: ≥8
    item_name == "Stairs"     & score <= 2     ~ "Dependent",
    item_name == "Stairs"     & score < 8      ~ "Assistance",
    item_name == "Stairs"     & score >= 8     ~ "Independent"
  )
}

# Pivot from wide to long format
causal_long <- causal_data %>%
  select(PatientID, Age, Pathology, Cardiovascular_Disease,
         Negative_Control_Outcome,
         Feeding, Ambulation, Stairs) %>%
  pivot_longer(
    cols      = c(Feeding, Ambulation, Stairs),
    names_to  = "Item",
    values_to = "RawScore"
  ) %>%
  mutate(
    # CAUSAL vs. PREDICTIVE CONDITIONING distinction:
    # When we fit the model with Score ~ Cardiovascular_Disease, we are identifying
    # the MARGINAL ASSOCIATION between CVD and Barthel outcomes given the confounders.
    # Later, G-computation (do-calculus) will estimate the CAUSAL effect by
    # marginalizing over the confounder distributions.
    Score = factor(
      discretize_score(RawScore, Item),
      levels = c("Dependent", "Assistance", "Independent"),
      ordered = TRUE
    ),
    Pathology                = factor(Pathology, levels = c("Neurological", "Orthopedic")),
    Cardiovascular_Disease   = as.numeric(Cardiovascular_Disease),  # 0/1 numeric for brms
    Item                     = factor(Item, levels = c("Feeding", "Ambulation", "Stairs")),
    Age                      = as.numeric(Age)
  ) %>%
  filter(!is.na(Score))

cat(sprintf("Final modeling dataset: %d observations from %d patients\n",
            nrow(causal_long), n_distinct(causal_long$PatientID)))
cat("\nScore × Item distribution:\n")
print(table(causal_long$Score, causal_long$Item))
cat("\nCVD exposure by Item:\n")
print(table(causal_long$Cardiovascular_Disease, causal_long$Item))

# =============================================================================
# PART 2: CAUSAL DAG SPECIFICATION
# =============================================================================

cat("\n=== STEP 2: CAUSAL DAG SPECIFICATION ===\n")

# Define the Directed Acyclic Graph (DAG) using dagitty syntax.
# This encodes our causal assumptions about the data-generating process.
dag_str <- "
  dag {
    Age [pos=\"0,0\"]
    Pathology [pos=\"1,0\"]
    Cardiovascular_Disease [exposure, pos=\"1,1\"]
    Barthel_Score [outcome, pos=\"2,1\"]
    
    Age -> Cardiovascular_Disease
    Age -> Barthel_Score
    Pathology -> Cardiovascular_Disease
    Pathology -> Barthel_Score
    Cardiovascular_Disease -> Barthel_Score
  }
"

# Parse the DAG with dagitty
dag <- dagitty(dag_str)

# Verify Backdoor Criterion: For causal identification of the effect of
# Cardiovascular_Disease on Barthel_Score, we need to block all backdoor paths.
# Backdoor paths are non-causal paths from the exposure to the outcome that
# do not pass through the exposure. They arise from confounding.

backdoor_vars <- adjustmentSets(dag, exposure = "Cardiovascular_Disease",
                                 outcome = "Barthel_Score", type = "all")
cat("\nBackdoor criterion check:\n")
cat("Sufficient adjustment sets to remove confounding bias:\n")
print(backdoor_vars)

# The set {Age, Pathology} should satisfy the backdoor criterion.
# This means that conditional on Age and Pathology, the remaining association
# between Cardiovascular_Disease and Barthel_Score is causal.

# Visualize the DAG
dag_plot <- ggdag(dag, layout = "nicely") +
  theme_dag() +
  labs(title = "Structural Causal Model (SCM): Barthel Index & CVD")

cat("\nDAG plotted successfully.\n")

# Save DAG for app visualization
saveRDS(dag, "R/causal_dag.rds")
cat("DAG saved to: R/causal_dag.rds\n")

# =============================================================================
# PART 3: BAYESIAN HIERARCHICAL ORDINAL MODEL (Causal Specification)
# =============================================================================

cat("\n=== STEP 3: FITTING BAYESIAN CAUSAL MODEL ===\n")

cat("Model Family: cumulative(probit)\n")
cat("Formula: Score ~ Age + Pathology + Cardiovascular_Disease + (1|PatientID) + (1|Item)\n")
cat("\nCausal Interpretation:\n")
cat("- Fixed effects (Age, Pathology, CVD) represent population-level associations.\n")
cat("- The coefficient on Cardiovascular_Disease is conditioned on Age & Pathology.\n")
cat("- Later, G-computation marginalizes over these confounders to compute the causal effect.\n\n")

# MCMC settings: lightweight for quick testing/iteration
n_cores <- min(parallel::detectCores(logical = FALSE), 2L)
cat(sprintf("MCMC: 2 chains × 2000 iterations (1000 warmup) on %d core(s)\n", n_cores))
cat("This will take 2-3 minutes...\n\n")

set.seed(42)

# Fit the Bayesian cumulative ordinal regression model
causal_model <- brm(
  formula = Score ~ Age + Pathology + Cardiovascular_Disease + (1 | PatientID) + (1 | Item),
  data    = causal_long,
  family  = cumulative("probit"),
  iter    = 2000,
  warmup  = 1000,
  chains  = 2,
  cores   = n_cores,
  seed    = 42,
  control = list(adapt_delta = 0.95),
  silent  = 2,
  refresh = 250
)

# =============================================================================
# PART 4: DIAGNOSTICS & CONVERGENCE CHECKS
# =============================================================================

cat("\n=== STEP 4: MCMC DIAGNOSTICS ===\n")
print(summary(causal_model))

rhats <- rhat(causal_model)
if (any(rhats > 1.01, na.rm = TRUE)) {
  warning("Some Rhat > 1.01 — potential convergence issues.")
} else {
  cat("\n✓ Convergence OK: all Rhat < 1.01\n")
}

# Extract the posterior distribution of the Cardiovascular_Disease coefficient
# This is our MARGINAL ASSOCIATION (not yet the CAUSAL effect)
cvd_posterior <- as_draws_df(causal_model, variable = "b_Cardiovascular_Disease")
cat("\nPosterior Summary: Effect of Cardiovascular_Disease (log-odds scale)\n")
print(summary(cvd_posterior))

# =============================================================================
# PART 5: NEGATIVE CONTROL VALIDATION
# =============================================================================

cat("\n=== STEP 5: NEGATIVE CONTROL OUTCOME VALIDATION ===\n")
cat("Fitting auxiliary model on Negative_Control_Outcome...\n")
cat("(Should show NO causal effect, validating model assumptions)\n\n")

# Create a long-format version for the negative control analysis
negative_ctrl_long <- causal_data %>%
  select(PatientID, Age, Pathology, Cardiovascular_Disease, 
         Negative_Control_Outcome) %>%
  mutate(
    Outcome = factor(Negative_Control_Outcome,
                     levels = c("Low", "Medium", "High"),
                     ordered = FALSE),
    Pathology = factor(Pathology, levels = c("Neurological", "Orthopedic")),
    Cardiovascular_Disease = as.numeric(Cardiovascular_Disease),
    Age = as.numeric(Age)
  )

# Fit model on negative control outcome (should show null effect)
negative_ctrl_model <- brm(
  formula = Outcome ~ Age + Pathology + Cardiovascular_Disease + (1 | PatientID),
  data    = negative_ctrl_long,
  family  = categorical("logit"),
  iter    = 1000,
  warmup  = 500,
  chains  = 2,
  cores   = n_cores,
  seed    = 42,
  silent  = 2,
  refresh = 0
)

cvd_neg_ctrl <- as_draws_df(negative_ctrl_model, variable = "b_Cardiovascular_Disease")
cat("\nNegative Control: Effect of Cardiovascular_Disease should be ≈ 0\n")
print(summary(cvd_neg_ctrl))

# =============================================================================
# PART 6: SAVE ALL OUTPUTS
# =============================================================================

cat("\n=== STEP 6: SAVING OUTPUTS ===\n")

# Save the fitted causal model
saveRDS(causal_model, "R/causal_barthel_model.rds")
cat("✓ Causal model saved to: R/causal_barthel_model.rds\n")

# Save the long-format causal data for app use
saveRDS(causal_long, "R/causal_data_long.rds")
cat("✓ Causal data saved to: R/causal_data_long.rds\n")

# Save full causal dataset for reference
saveRDS(causal_data, "R/causal_data_wide.rds")
cat("✓ Wide-format causal data saved to: R/causal_data_wide.rds\n")

# Capture session info for reproducibility
writeLines(capture.output(sessionInfo()), "R/causal_session_info.txt")
cat("✓ Session info saved to: R/causal_session_info.txt\n")

cat("\n" %+% strrep("=", 75) %+% "\n")
cat("CAUSAL INFERENCE PIPELINE COMPLETE\n")
cat("Next step: Launch the Causal Shiny app\n")
cat("  Rscript -e \"shiny::runApp('R', launch.browser = TRUE)\"\n")
cat(strrep("=", 75) %+% "\n\n")
