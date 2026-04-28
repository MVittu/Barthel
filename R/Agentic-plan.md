You are a Senior Biostatistician and an Advanced R/Shiny Developer. Your task is to create a complete infrastructure to analyze the "Seccia et al. 2020" dataset (clinical data and Barthel indices) and produce an interactive Clinical Decision Support Dashboard in Shiny.

The scientific objective is to overcome the "metric fallacy" of the Barthel Index total score by modeling the individual items (Feeding, Ambulation, Stairs, etc.) using a Bayesian Hierarchical Ordinal Model. This will provide individualized recovery probability estimates with 95% Credible Intervals (95% CrI).

You must generate two distinct R scripts:
1. `01_model_training.R`: For preprocessing and fitting the Bayesian model.
2. `app.R`: For the interactive Shiny user interface.

Technical Requirements and Libraries:
- Data Wrangling: `dplyr`, `tidyr`
- Bayesian Modeling: `brms`
- Uncertainty Extraction: `tidybayes`
- Visualization: `ggplot2`, `bslib` (for a modern UI theme)

=========================================
TASK 1: Creation of `01_model_training.R`
=========================================
1. Simulate a mock dataset representative of "Seccia et al. 2020" (if you do not have access to the real file). The dataset must be in wide format with the following columns: PatientID, Age, Pathology (Neurological/Orthopedic), Comorbidity (e.g., None, Cardiovascular), and the discharge Barthel columns (COD_26_1 for Feeding, COD_26_9 for Ambulation, COD_26_10 for Stairs). The Barthel values must be ordinal (e.g., 0, 5, 10, 15).
2. Execute a `pivot_longer` to transform the Barthel columns into Long format (columns: PatientID, Age, Pathology, Comorbidity, Item, Score). Transform "Score" into an ordered factor.
3. Configure and fit a Bayesian ordinal regression model using `brms`. 
   - Suggested formula: `Score ~ Age + Pathology + Comorbidity + (1 | PatientID) + (1 | Item)`
   - Family: `cumulative("probit")` or `cumulative("logit")`.
   - Set lightweight computation parameters (e.g., iter = 2000, chains = 2) so the code can be tested quickly.
4. Save the fitted model as an `.rds` file (e.g., `barthel_model.rds`) in the current directory.

=========================================
TASK 2: Creation of `app.R` (Shiny App)
=========================================
The app must load `barthel_model.rds` at startup and use it to perform dynamic inference.

UI (User Interface) using `bslib`:
- Layout: Sidebar layout.
- Sidebar: User-controlled inputs.
  - `selectInput` for "Main Pathology" (Neurological, Orthopedic).
  - `sliderInput` for "Age" (from 30 to 90 years).
  - `selectInput` for "Comorbidity" (None, Cardiovascular, Diabetes).
- Main Panel: 
  - Three side-by-side plots (or stacked in a column if space is insufficient) dedicated to: Feeding, Ambulation, Stairs.
  - Below the plots, a `uiOutput` or text block (bslib Card) for the "Biostatistical Synthesis and Recommendations".

Server Logic:
1. Reactivity: Create a `reactive` dataframe containing a "new patient" based on the current UI inputs, expanded for the three target Items (Feeding, Ambulation, Stairs).
2. Bayesian Inference: Use `add_epred_draws()` from `tidybayes` combined with `fitted()` or `posterior_epred()` from `brms` to generate the posterior probability distribution for each score category.
3. Uncertainty Calculation: Group by Item and Score Category, calculating the mean (Expected Probability) and the lower and upper bounds at 95% CrI (e.g., 0.025 and 0.975 quantiles).
4. Visualization (`ggplot2`): 
   - Generate vertical bar charts for each Item. The X-axis is the score category (Dependent, Assistance, Independent), and the Y-axis is the probability (0-100%).
   - YOU MUST ADD `geom_errorbar` to show the 95% CrI calculated in step 3.
   - Ensure the plots have a clean theme (`theme_minimal()`) and clear titles.
5. Dynamic Synthesis: Write a reactive logic that generates a brief interpretive text. For example, if the credible interval for independence in ambulation is very wide (difference > 30%), the text should print a warning about high prognostic uncertainty.

Required Output: 
Generate complete, clean, and extensively commented R code for both files. Handle any potential errors and provide clear instructions on how to run everything in VS Code.