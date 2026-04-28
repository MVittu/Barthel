# =============================================================================
# app.R
# Barthel Index — Clinical Decision Support Dashboard
# Bayesian Hierarchical Ordinal Model (brms) — Interactive Shiny Application
#
# Prerequisites:
#   1. Run 01_model_training.R to generate barthel_model.rds
#   2. Required packages: shiny, bslib, dplyr, ggplot2, brms, tidybayes
#
# How to run in VS Code terminal:
#   Rscript -e "shiny::runApp('.')"
# or in R/RStudio:
#   shiny::runApp()
# =============================================================================

# ── 0. LIBRARIES ──────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(dplyr)
  library(ggplot2)
  library(brms)
  library(tidybayes)
})

# ── 1. LOAD MODEL AT STARTUP (shared across all sessions) ─────────────────────
if (!file.exists("barthel_model.rds")) {
  stop(
    "barthel_model.rds not found.\n",
    "Please run 01_model_training.R first to fit and save the model."
  )
}

cat("Loading Bayesian model… ")
bart_model <- readRDS("barthel_model.rds")
cat("done.\n")

# ── 2. USER INTERFACE ─────────────────────────────────────────────────────────
ui <- page_sidebar(
  title = tags$span(
    icon("heartbeat"), " Barthel Index · Clinical Decision Support Dashboard"
  ),
  theme = bs_theme(
    bootswatch = "flatly",
    primary    = "#2C3E50",
    secondary  = "#18BC9C"
  ),

  # ── Sidebar: patient profile inputs ────────────────────────────────────────
  sidebar = sidebar(
    width = 280,
    h5(icon("user-md"), " Patient Profile", class = "text-primary"),
    hr(),

    selectInput(
      inputId  = "pathology",
      label    = "Main Pathology",
      choices  = c("Neurological", "Orthopedic"),
      selected = "Neurological"
    ),
    sliderInput(
      inputId = "age",
      label   = "Age (years)",
      min     = 30, max = 90, value = 65, step = 1
    ),
    selectInput(
      inputId  = "comorbidity",
      label    = "Comorbidity",
      choices  = c("None", "Cardiovascular", "Diabetes"),
      selected = "None"
    ),
    hr(),
    actionButton(
      inputId = "predict_btn",
      label   = tagList(icon("sync"), " Update Predictions"),
      class   = "btn-primary w-100"
    ),
    hr(),
    p(class = "small text-muted",
      "Bayesian Hierarchical Ordinal Model (brms, cumulative-probit).",
      "Bars show posterior mean probability; error bars show 95% Credible Intervals.",
      "Initial predictions load automatically.")
  ),

  # ── Main panel: three item plots ────────────────────────────────────────────
  layout_columns(
    col_widths = c(4, 4, 4),
    card(
      card_header(icon("utensils"), " Feeding"),
      card_body(plotOutput("plot_feeding", height = "290px"))
    ),
    card(
      card_header(icon("walking"), " Ambulation"),
      card_body(plotOutput("plot_ambulation", height = "290px"))
    ),
    card(
      card_header(icon("stairs"), " Stairs"),
      card_body(plotOutput("plot_stairs", height = "290px"))
    )
  ),

  # ── Biostatistical synthesis card ──────────────────────────────────────────
  card(
    card_header(
      class = "bg-primary text-white",
      icon("stethoscope"), " Biostatistical Synthesis & Clinical Recommendations"
    ),
    card_body(uiOutput("synthesis"))
  )
)

# ── 3. SERVER LOGIC ───────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── 3.1 Reactive: new-patient data frame ────────────────────────────────────
  # Depends on predict_btn; isolate() prevents re-running on each slider move.
  # Fires on startup (button value = 0) and on each click.
  new_pt <- reactive({
    input$predict_btn  # reactive dependency — re-evaluate when button is clicked
    isolate({
      tibble(
        PatientID   = "NEW_PATIENT",
        Age         = as.numeric(input$age),
        # Factor levels must match those used during model training
        Pathology   = factor(input$pathology,
                             levels = c("Neurological", "Orthopedic")),
        Comorbidity = factor(input$comorbidity,
                             levels = c("Cardiovascular", "Diabetes", "None")),
        Item        = factor(c("Feeding", "Ambulation", "Stairs"),
                             levels = c("Feeding", "Ambulation", "Stairs"))
      )
    })
  })

  # ── 3.2 Reactive: posterior predictive summary ──────────────────────────────
  # add_epred_draws() returns draws of the expected (mean) category probabilities.
  # allow_new_levels = TRUE lets brms handle "NEW_PATIENT" not seen in training
  # by sampling its random effect from the population prior.
  post_summary <- reactive({
    pt <- new_pt()

    withProgress(message = "Computing Bayesian predictions…", value = 0, {
      setProgress(0.3)

      draws <- tryCatch(
        pt %>%
          add_epred_draws(
            object           = bart_model,
            allow_new_levels = TRUE,
            ndraws           = 500,   # 500 draws → fast but smooth CrI
            seed             = 42,
            re_formula       = NULL   # include all random effects (patient + item)
          ),
        error = function(e) {
          showNotification(paste("Prediction error:", conditionMessage(e)),
                           type = "error", duration = 8)
          return(NULL)
        }
      )
      setProgress(1)
    })

    req(draws)

    draws %>%
      group_by(Item, .category) %>%
      summarise(
        mean_p = mean(.epred),
        lo95   = quantile(.epred, 0.025),
        hi95   = quantile(.epred, 0.975),
        .groups = "drop"
      ) %>%
      mutate(
        .category = factor(.category,
                           levels = c("Dependent", "Assistance", "Independent"))
      )
  })

  # ── 3.3 Helper: build a single-item probability bar chart ───────────────────
  #   X-axis : recovery level (Dependent / Assistance / Independent)
  #   Y-axis : posterior mean probability (0–100%)
  #   Error bars: 95% Credible Interval
  make_plot <- function(item_name, fill_hex) {
    d <- post_summary() %>% filter(Item == item_name)

    ggplot(d, aes(x = .category, y = mean_p * 100)) +
      geom_col(fill = fill_hex, alpha = 0.85, width = 0.55) +
      geom_errorbar(
        aes(ymin = lo95 * 100, ymax = hi95 * 100),
        width     = 0.22,
        linewidth = 0.9,
        colour    = "grey30"
      ) +
      geom_text(
        aes(label = sprintf("%.1f%%", mean_p * 100)),
        vjust      = -0.7,
        size       = 3.6,
        fontface   = "bold"
      ) +
      scale_y_continuous(
        limits = c(0, 110),
        labels = function(x) paste0(x, "%"),
        expand = expansion(mult = c(0, 0))
      ) +
      labs(x = "Recovery Level", y = "Probability (95% CrI)", title = item_name) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title         = element_text(hjust = 0.5, face = "bold", size = 13),
        axis.title         = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank()
      )
  }

  # ── 3.4 Render plots ────────────────────────────────────────────────────────
  output$plot_feeding <- renderPlot({
    req(post_summary())
    make_plot("Feeding",    "#2980B9")
  })
  output$plot_ambulation <- renderPlot({
    req(post_summary())
    make_plot("Ambulation", "#27AE60")
  })
  output$plot_stairs <- renderPlot({
    req(post_summary())
    make_plot("Stairs",     "#8E44AD")
  })

  # ── 3.5 Dynamic synthesis text ──────────────────────────────────────────────
  # Generates interpretive text with clinical alerts, using CrI width as a proxy
  # for prognostic uncertainty.
  output$synthesis <- renderUI({
    req(post_summary())
    s <- post_summary()

    # Extract independence statistics for each item
    item_stats <- function(item) {
      r <- s %>% filter(Item == item, .category == "Independent")
      list(m = r$mean_p, lo = r$lo95, hi = r$hi95,
           w = r$hi95 - r$lo95)  # CrI width (proxy for uncertainty)
    }
    f  <- item_stats("Feeding")
    am <- item_stats("Ambulation")
    st <- item_stats("Stairs")

    mean_indep <- mean(c(f$m, am$m, st$m))

    # Overall prognosis tier
    overall <- if (mean_indep >= 0.60) {
      list(txt = "Favourable Overall Prognosis",  cls = "success")
    } else if (mean_indep >= 0.35) {
      list(txt = "Moderate Overall Prognosis",    cls = "warning")
    } else {
      list(txt = "Guarded Overall Prognosis",     cls = "danger")
    }

    # Uncertainty badge per item
    precision_badge <- function(width) {
      if (width > 0.30)
        tags$span(class = "badge bg-warning text-dark", "High Uncertainty")
      else if (width > 0.15)
        tags$span(class = "badge bg-info text-white",  "Moderate Uncertainty")
      else
        tags$span(class = "badge bg-success",           "Good Precision")
    }

    # Clinical recommendation bullets
    notes <- list()

    # Flag: wide CrI on Ambulation → recommend re-evaluation
    if (am$w > 0.30)
      notes <- c(notes, list(tags$li(
        tags$strong("Ambulation:"),
        sprintf(
          " Wide 95%% CrI (%.0f percentage points). Prognostic uncertainty is high.",
          am$w * 100),
        " Recommend re-evaluation after 2 weeks of rehabilitation."
      )))

    # Flag: low independence probability for Stairs → consider home adaptation
    if (st$m < 0.25)
      notes <- c(notes, list(tags$li(
        tags$strong("Stairs:"),
        sprintf(
          " Low probability of independent stair use (%.1f%%).", st$m * 100),
        " Consider targeted stair-training programme or home adaptation assessment."
      )))

    # Flag: high independence probability for Feeding → reassurance
    if (f$m > 0.75)
      notes <- c(notes, list(tags$li(
        tags$strong("Feeding:"),
        sprintf(
          " High probability of independent feeding (%.1f%%).", f$m * 100),
        " Standard nutritional supervision should suffice."
      )))

    # Flag: wide CrI on Feeding
    if (f$w > 0.30)
      notes <- c(notes, list(tags$li(
        tags$strong("Feeding:"),
        sprintf(
          " Considerable uncertainty in feeding outcome (CrI width %.0f pp).",
          f$w * 100),
        " Regular occupational therapy assessment recommended."
      )))

    if (length(notes) == 0)
      notes <- list(tags$li("No specific clinical alerts for this patient profile."))

    # Compose HTML output
    tagList(
      # Overall alert banner
      div(
        class = paste0("alert alert-", overall$cls, " mb-3"),
        tags$strong(overall$txt),
        sprintf(" — Mean P(Independence) across items: %.1f%%", mean_indep * 100)
      ),

      # Summary table: one row per Barthel item
      tags$table(
        class = "table table-sm table-bordered table-striped mb-3",
        tags$thead(
          tags$tr(
            tags$th("Barthel Item"),
            tags$th("P(Independent)"),
            tags$th("95% Credible Interval"),
            tags$th("Precision")
          )
        ),
        tags$tbody(
          lapply(
            list(list("Feeding", f), list("Ambulation", am), list("Stairs", st)),
            function(x) {
              nm <- x[[1]]; v <- x[[2]]
              tags$tr(
                tags$td(tags$strong(nm)),
                tags$td(sprintf("%.1f%%", v$m  * 100)),
                tags$td(sprintf("[%.1f%% – %.1f%%]", v$lo * 100, v$hi * 100)),
                tags$td(precision_badge(v$w))
              )
            }
          )
        )
      ),

      # Clinical notes
      tags$div(
        tags$strong("Clinical Notes:"),
        tags$ul(notes),
        tags$p(
          class = "small text-muted mt-2",
          sprintf(
            "Patient profile: %s · Age %d · Comorbidity: %s.",
            isolate(input$pathology), isolate(input$age), isolate(input$comorbidity)
          )
        )
      )
    )
  })
}

# ── 4. LAUNCH APPLICATION ─────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
