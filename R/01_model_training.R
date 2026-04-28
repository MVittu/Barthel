# =============================================================================
# 01_model_training.R
# Barthel Index — Bayesian Hierarchical Ordinal Model Training
# Reference: Seccia et al. 2020 Dataset
#
# Scientific goal: overcome the "metric fallacy" of the Barthel total score by
# modelling individual items (Feeding, Ambulation, Stairs) with a Bayesian
# Hierarchical Ordinal Model, producing recovery-probability estimates with
# 95% Credible Intervals.
#
# Output: barthel_model.rds  (loaded by app.R)
#
# How to run in VS Code terminal:
#   Rscript 01_model_training.R
# or interactively:
#   source("01_model_training.R")
# =============================================================================

# ── 0. LIBRARIES ──────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(brms)
  library(tidybayes)
})

# ── 1. ICD-9 CLASSIFICATION HELPERS ──────────────────────────────────────────

# Classifies primary ICD-9-CM code into "Neurological", "Orthopedic", or "Other".
# Matching logic is prefix-based (codes stored without decimal points, e.g.
# ICD 438.21 → "43821").
classify_pathology <- function(code) {
  s <- trimws(toupper(as.character(code)))
  if (s %in% c("99", "999", "", "NA")) return("Other")

  # Orthopedic: joint replacement V43.6x, post-arthroplasty V45.x
  if (grepl("^V43[6-9]|^V450|^V451", s)) return("Orthopedic")

  # Neurological: cerebrovascular disease (430–438)
  if (grepl("^43[0-8]", s)) return("Neurological")
  # Neurological: CNS & peripheral nervous system (320–359)
  if (grepl("^3[2-5][0-9]", s)) return("Neurological")
  # Neurological: neurological signs/symptoms incl. gait disorders (780–781)
  if (grepl("^78[01]", s)) return("Neurological")

  # Orthopedic: musculoskeletal & connective tissue (710–739)
  if (grepl("^7[1-3][0-9]", s)) return("Orthopedic")
  # Orthopedic: fractures of skull, spine, pelvis, limbs (800–829)
  if (grepl("^8[0-2][0-9]", s)) return("Orthopedic")

  return("Other")
}

# Scans a list of secondary ICD-9 codes for cardiovascular or diabetes
# comorbidities.  Returns "Cardiovascular", "Diabetes", or "None".
classify_comorbidity <- function(code_list) {
  has_cv <- FALSE
  has_dm <- FALSE
  for (code in code_list) {
    s <- trimws(toupper(as.character(code)))
    if (is.na(s) || s %in% c("99", "999", "")) next
    # Diabetes mellitus (250–259)
    if (grepl("^25[0-9]", s))                      has_dm <- TRUE
    # Cardiovascular: hypertension 401–405, IHD 410–414,
    # other heart 420–429, arteries 440–448, AF (427)
    if (grepl("^40[1-5]|^41[0-4]|^42[0-9]|^44[0-8]|^427", s)) has_cv <- TRUE
  }
  if (has_cv) return("Cardiovascular")
  if (has_dm) return("Diabetes")
  "None"
}

# Maps a raw Barthel item score to a 3-level ordered category.
# Thresholds follow established clinical cut-points for each item:
#   Feeding    (max 10): 0 = Dependent, 5 = Assistance, ≥8 = Independent
#   Ambulation (max 15): ≤3 = Dependent, ≤11 = Assistance, ≥12 = Independent
#   Stairs     (max 10): ≤2 = Dependent, ≤7 = Assistance, ≥8 = Independent
categorize_score <- function(score, item) {
  x <- suppressWarnings(as.numeric(score))
  if (is.na(x)) return(NA_character_)
  switch(item,
    Feeding    = if (x == 0) "Dependent" else if (x < 8)  "Assistance" else "Independent",
    Ambulation = if (x <= 3) "Dependent" else if (x < 12) "Assistance" else "Independent",
    Stairs     = if (x <= 2) "Dependent" else if (x < 8)  "Assistance" else "Independent",
    NA_character_
  )
}

# ── 2. DATA LOADING ───────────────────────────────────────────────────────────
cat("=== Step 1: Loading dataset ===\n")
raw <- read.csv("data/Dataset.csv",
                stringsAsFactors = FALSE,
                na.strings       = c("", "NA"))

# Strip leading/trailing whitespace introduced by spaces in the CSV header
names(raw) <- trimws(names(raw))

cat(sprintf("Loaded: %d rows × %d columns\n", nrow(raw), ncol(raw)))

# ── 3. FEATURE ENGINEERING ────────────────────────────────────────────────────
cat("\n=== Step 2: Classifying Pathology & Comorbidity ===\n")

# Secondary diagnosis columns (COD_2 … COD_10) used for comorbidity search
sec_diag_cols <- paste0("COD_", 2:10)

proc <- raw %>%
  mutate(
    PatientID = paste0("PT", row_number()),
    # Primary diagnosis → Pathology group
    Pathology = sapply(COD_1, classify_pathology),
    # Secondary diagnoses → dominant comorbidity
    Comorbidity = apply(
      .[, sec_diag_cols], 1,
      function(r) classify_comorbidity(as.list(r))
    )
  ) %>%
  # Keep only rehabilitation-relevant pathology groups
  filter(Pathology %in% c("Neurological", "Orthopedic"))

cat("Pathology distribution (after filtering):\n")
print(table(proc$Pathology))
cat("\nComorbidity distribution:\n")
print(table(proc$Comorbidity))

# ── 4. BARTHEL ITEMS — WIDE → LONG FORMAT ────────────────────────────────────
cat("\n=== Step 3: Pivot Barthel items to long format ===\n")

# Target discharge Barthel items:
#   COD_26_1  = Feeding
#   COD_26_9  = Ambulation
#   COD_26_10 = Stairs
long <- proc %>%
  select(PatientID, Age, Pathology, Comorbidity,
         Feeding    = COD_26_1,
         Ambulation = COD_26_9,
         Stairs     = COD_26_10) %>%
  pivot_longer(
    cols      = c(Feeding, Ambulation, Stairs),
    names_to  = "Item",
    values_to = "RawScore"
  ) %>%
  mutate(
    # Convert raw numeric score → ordered 3-level factor
    Score = factor(
      mapply(categorize_score, RawScore, Item),
      levels  = c("Dependent", "Assistance", "Independent"),
      ordered = TRUE
    ),
    # Set factor levels explicitly so all levels appear in the model
    Pathology   = factor(Pathology,   levels = c("Neurological", "Orthopedic")),
    Comorbidity = factor(Comorbidity, levels = c("Cardiovascular", "Diabetes", "None")),
    Item        = factor(Item,        levels = c("Feeding", "Ambulation", "Stairs")),
    Age         = as.numeric(Age)
  ) %>%
  filter(!is.na(Score), !is.na(Age))

cat(sprintf("Final modelling dataset: %d observations from %d patients\n",
            nrow(long), n_distinct(long$PatientID)))
cat("\nScore × Item distribution:\n")
print(table(long$Score, long$Item))

# ── 5. BAYESIAN ORDINAL MODEL ─────────────────────────────────────────────────
cat("\n=== Step 4: Fitting Bayesian Hierarchical Ordinal Model ===\n")
cat("Family  : cumulative(probit)\n")
cat("Formula : Score ~ Age + Pathology + Comorbidity + (1|PatientID) + (1|Item)\n")
cat("Chains  : 2   |  Iterations: 2000  |  Warmup: 1000\n")
cat("This will take several minutes — please wait...\n\n")

# cumulative("probit") models the probability of being in each ordered category.
# Random intercepts per patient and per item capture between-patient and
# between-item baseline variability (the hierarchical structure).
bart_model <- brm(
  formula = Score ~ Age + Pathology + Comorbidity + (1 | PatientID) + (1 | Item),
  data    = long,
  family  = cumulative("probit"),
  iter    = 2000,
  warmup  = 1000,
  chains  = 2,
  cores   = 2,
  seed    = 42,
  control = list(adapt_delta = 0.95),
  silent  = 2,
  refresh = 250
)

# ── 6. DIAGNOSTICS ────────────────────────────────────────────────────────────
cat("\n=== Step 5: Diagnostics ===\n")
cat("\nModel Summary:\n")
print(summary(bart_model))

rhats <- rhat(bart_model)
if (any(rhats > 1.01, na.rm = TRUE)) {
  warning("Some Rhat > 1.01 — potential convergence issues. Consider increasing iter/chains.")
} else {
  cat("\nConvergence OK: all Rhat < 1.01\n")
}

# ── 7. SAVE ───────────────────────────────────────────────────────────────────
cat("\n=== Step 6: Saving model ===\n")
saveRDS(bart_model, "barthel_model.rds")
cat("Model saved to: barthel_model.rds\n")
cat("Next step    : open app.R and run shiny::runApp() or press Run App in RStudio.\n")
