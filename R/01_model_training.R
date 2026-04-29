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
# How to run in VS Code terminal (from project root):
#   Rscript R/01_model_training.R
# or interactively (working directory = project root):
#   source("R/01_model_training.R")
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(brms)
  library(tidybayes)
})

# ── WORKING DIRECTORY GUARD ───────────────────────────────────────────────────
# All relative paths assume the project root as CWD. If the script is invoked
# from inside R/ (e.g. cd R && Rscript 01_model_training.R), step up one level.
if (!file.exists("data/Dataset.csv")) {
  if (file.exists(file.path("..", "data", "Dataset.csv"))) {
    setwd("..")
  } else {
    stop(
      "Cannot locate data/Dataset.csv.\n",
      "Run from the project root: Rscript R/01_model_training.R\n",
      "Current directory: ", normalizePath(getwd())
    )
  }
}

# ── ICD-9 CLASSIFICATION HELPERS ─────────────────────────────────────────────

normalize_icd <- function(code) trimws(toupper(as.character(code)))

# Prefix-based: codes stored without decimal points (e.g. ICD 438.21 → "43821").
classify_pathology <- function(code) {
  s <- normalize_icd(code)
  if (s %in% c("99", "999", "", "NA")) return("Other")
  # Orthopedic: joint replacement V43.6x, post-arthroplasty V45.x
  if (grepl("^V43[6-9]|^V450|^V451", s)) return("Orthopedic")
  # Neurological: cerebrovascular (430–438), CNS/PNS (320–359), symptoms (780–781)
  if (grepl("^43[0-8]|^3[2-5][0-9]|^78[01]", s)) return("Neurological")
  # Orthopedic: musculoskeletal (710–739), fractures (800–829)
  if (grepl("^7[1-3][0-9]|^8[0-2][0-9]", s)) return("Orthopedic")
  "Other"
}

# Scans secondary ICD-9 codes for cardiovascular (401–448) or diabetes (250–259).
classify_comorbidity <- function(code_list) {
  has_cv <- FALSE
  has_dm <- FALSE
  for (code in code_list) {
    s <- normalize_icd(code)
    if (is.na(s) || s %in% c("99", "999", "")) next
    if (!has_dm) has_dm <- grepl("^25[0-9]", s)
    if (!has_cv) has_cv <- grepl("^40[1-5]|^41[0-4]|^42[0-9]|^44[0-8]|^427", s)
    if (has_cv && has_dm) break
  }
  if (has_cv) return("Cardiovascular")
  if (has_dm) return("Diabetes")
  "None"
}

# ── 1. DATA LOADING ───────────────────────────────────────────────────────────
cat("=== Step 1: Loading dataset ===\n")
raw <- read.csv("data/Dataset.csv",
                stringsAsFactors = FALSE,
                na.strings       = c("", "NA"))
names(raw) <- trimws(names(raw))
cat(sprintf("Loaded: %d rows × %d columns\n", nrow(raw), ncol(raw)))

# ── 2. FEATURE ENGINEERING ────────────────────────────────────────────────────
cat("\n=== Step 2: Classifying Pathology & Comorbidity ===\n")

sec_diag_cols <- paste0("COD_", 2:10)

proc <- raw %>%
  mutate(
    PatientID   = paste0("PT", row_number()),
    Pathology   = sapply(COD_1, classify_pathology),
    Comorbidity = apply(.[, sec_diag_cols], 1,
                        function(r) classify_comorbidity(as.list(r)))
  ) %>%
  filter(Pathology %in% c("Neurological", "Orthopedic"))

cat("Pathology distribution (after filtering):\n")
print(table(proc$Pathology))
cat("\nComorbidity distribution:\n")
print(table(proc$Comorbidity))

# ── 3. BARTHEL ITEMS — WIDE → LONG FORMAT ────────────────────────────────────
cat("\n=== Step 3: Pivot Barthel items to long format ===\n")

# Discharge Barthel items: COD_26_1 = Feeding, COD_26_9 = Ambulation, COD_26_10 = Stairs
long <- proc %>%
  select(PatientID, Age, Pathology, Comorbidity,
         Feeding    = COD_26_1,
         Ambulation = COD_26_9,
         Stairs     = COD_26_10) %>%
  pivot_longer(c(Feeding, Ambulation, Stairs), names_to = "Item", values_to = "RawScore") %>%
  mutate(
    RawScore = suppressWarnings(as.numeric(RawScore)),
    # Clinical cut-points per item:
    #   Feeding (max 10):    0=Dependent  | <8=Assistance  | ≥8=Independent
    #   Ambulation (max 15): ≤3=Dependent | <12=Assistance | ≥12=Independent
    #   Stairs (max 10):     ≤2=Dependent | <8=Assistance  | ≥8=Independent
    Score = factor(
      case_when(
        is.na(RawScore)                          ~ NA_character_,
        Item == "Feeding"    & RawScore == 0     ~ "Dependent",
        Item == "Feeding"    & RawScore < 8      ~ "Assistance",
        Item == "Feeding"                         ~ "Independent",
        Item == "Ambulation" & RawScore <= 3     ~ "Dependent",
        Item == "Ambulation" & RawScore < 12     ~ "Assistance",
        Item == "Ambulation"                      ~ "Independent",
        Item == "Stairs"     & RawScore <= 2     ~ "Dependent",
        Item == "Stairs"     & RawScore < 8      ~ "Assistance",
        Item == "Stairs"                          ~ "Independent"
      ),
      levels = c("Dependent", "Assistance", "Independent"),
      ordered = TRUE
    ),
    # Explicit levels ensure brms sees all categories even if sparse in data
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

# ── 4. BAYESIAN ORDINAL MODEL ─────────────────────────────────────────────────
cat("\n=== Step 4: Fitting Bayesian Hierarchical Ordinal Model ===\n")
cat("Family  : cumulative(probit)\n")
cat("Formula : Score ~ Age + Pathology + Comorbidity + (1|PatientID) + (1|Item)\n")

# Detect available physical cores; cap at 2 so the script runs on any machine.
n_cores <- min(parallel::detectCores(logical = FALSE), 2L)
cat(sprintf("Chains  : 2   |  Iterations: 2000  |  Warmup: 1000  |  Cores: %d\n", n_cores))
cat("This will take several minutes — please wait...\n\n")

# set.seed covers any R-level randomness before Stan takes over with seed = 42.
set.seed(42)

bart_model <- brm(
  formula = Score ~ Age + Pathology + Comorbidity + (1 | PatientID) + (1 | Item),
  data    = long,
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

# ── 5. DIAGNOSTICS ────────────────────────────────────────────────────────────
cat("\n=== Step 5: Diagnostics ===\n")
print(summary(bart_model))

rhats <- rhat(bart_model)
if (any(rhats > 1.01, na.rm = TRUE)) {
  warning("Some Rhat > 1.01 — potential convergence issues. Consider increasing iter/chains.")
} else {
  cat("\nConvergence OK: all Rhat < 1.01\n")
}

# ── 6. SAVE ───────────────────────────────────────────────────────────────────
cat("\n=== Step 6: Saving model & session info ===\n")
saveRDS(bart_model, "barthel_model.rds")
cat("Model saved to: barthel_model.rds\n")

# Session snapshot lets any collaborator reproduce the exact environment.
# Note: Stan/MCMC results are statistically equivalent across platforms with
# seed = 42 but may not be bit-for-bit identical due to floating-point
# differences between OS/compiler/BLAS implementations.
writeLines(capture.output(sessionInfo()), "session_info.txt")
cat("Session info saved to: session_info.txt\n")
cat("Next step: Rscript -e \"shiny::runApp('R')\"\n")
