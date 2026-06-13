library(shiny)
library(ggplot2)
library(dplyr)
library(tibble)
library(airqmra)
if (!requireNamespace("tidyr", quietly = TRUE)) {
  stop("Package 'tidyr' is needed for this app. Please install it with install.packages(\"tidyr\").")
}

# Pathogen list
pathogen_choices <- unique(airqmra::pathogen_db$pathogen_id)
names(pathogen_choices) <- paste0(
  airqmra::pathogen_db$pathogen_cn[match(pathogen_choices, airqmra::pathogen_db$pathogen_id)],
  " (", pathogen_choices, ")"
)

# UI
ui <- navbarPage(
  title = "airqmra Inhalation Risk Assessment",
  theme = shinythemes::shinytheme("flatly"),

  tabPanel("Scenario Settings",
           sidebarLayout(
             sidebarPanel(
               selectInput("pathogen", "Pathogen", choices = pathogen_choices),
               selectInput("exposure_group", "Population Group", choices = NULL),
               selectInput("microenvironment", "Microenvironment",
                           choices = c("Indoor", "Outdoor")),
               selectInput("season", "Season",
                           choices = c("Annual", "Spring & Autumn", "Summer", "Winter")),
               selectInput("population_deposition", "Deposition Parameter Group",
                           choices = c("Male", "Female", "Child", "Senior")),
               radioButtons("use_uncertainty", "Dose-Response Parameter Uncertainty",
                            choices = c("Point Estimate" = FALSE, "Include Uncertainty" = TRUE),
                            selected = FALSE),
               hr(),
               h5("Air Concentration (copies/m^3)"),
               radioButtons("conc_type", "Concentration Type",
                            choices = c("Default Lognormal" = "default",
                                        "Fixed Value" = "fixed",
                                        "Upload File" = "file"),
                            selected = "default"),
               conditionalPanel(
                 condition = "input.conc_type == 'fixed'",
                 numericInput("conc_fixed", "Fixed Concentration", value = 50, min = 0)
               ),
               conditionalPanel(
                 condition = "input.conc_type == 'file'",
                 fileInput("conc_file", "Choose CSV File",
                           accept = c("text/csv", ".csv")),
                 selectInput("conc_col", "Concentration Column", choices = NULL),
                 radioButtons("conc_method", "Fitting Method",
                              choices = c("Empirical" = "empirical",
                                          "Lognormal" = "lognormal",
                                          "Auto Select" = "auto"),
                              selected = "lognormal")
               ),
               hr(),
               h5("Annual Exposure Parameters"),
               radioButtons("annual_type", "Annual Frequency",
                            choices = c("Fixed Count" = "fixed_freq", "Poisson Distribution" = "poisson"),
                            selected = "fixed_freq"),
               conditionalPanel(
                 condition = "input.annual_type == 'fixed_freq'",
                 numericInput("freq", "Exposures per Year", value = 20, min = 1)
               ),
               conditionalPanel(
                 condition = "input.annual_type == 'poisson'",
                 numericInput("lambda", "Poisson Mean", value = 20, min = 0)
               ),
               hr(),
               h5("Simulation Settings"),
               numericInput("n_iter", "1D Simulation Iterations", value = 5000, min = 1000, step = 1000),
               numericInput("n_var", "2D - Variability Iterations (n_var)", value = 1000, min = 500, step = 500),
               numericInput("n_unc", "2D - Uncertainty Iterations (n_unc)", value = 300, min = 100, step = 100),
               hr(),
               actionButton("run1d", "Run 1D Simulation", icon = icon("play"), class = "btn-primary"),
               br(), br(),
               actionButton("run2d", "Run 2D Simulation", icon = icon("random"), class = "btn-success"),
               br(), br(),
               p("2D simulation is time-consuming. Please be patient."),
               hr(),
               h5("Report Generation"),
               radioButtons("report_format", "Format",
                            choices = c("HTML" = "html", "Word" = "word"), selected = "html"),
               radioButtons("report_lang", "Language",
                            choices = c("Chinese" = "zh", "English" = "en"), selected = "en"),
               textInput("report_dir", "Output Directory", value = getwd()),
               actionButton("set_dir", "Set to Current R Working Directory"),
               br(),
               actionButton("generate_report", "Generate Report", icon = icon("file-alt"))
             ),
             mainPanel(
               h4("Current Scenario Summary"),
               verbatimTextOutput("scenario_summary"),
               h4("Scenario Validation"),
               verbatimTextOutput("scenario_check")
             )
           )
  ),

  tabPanel("1D Simulation",
           plotOutput("hist1d"),
           verbatimTextOutput("summary1d"),
           h4("Annual Risk"),
           plotOutput("annual_hist"),
           verbatimTextOutput("annual_summary")
  ),

  tabPanel("2D Simulation",
           plotOutput("hist2d_overall"),
           verbatimTextOutput("summary2d"),
           h4("Uncertainty in Mean Risk"),
           plotOutput("hist2d_mean"),
           h4("Variability Distribution across Uncertainty Iterations"),
           plotOutput("var_boxplot")
  ),

  tabPanel("Sensitivity Analysis",
           h4("Spearman Tornado Plot (1D)"),
           plotOutput("tornado1d"),
           h4("Parameter Uncertainty (requires 2D simulation)"),
           actionButton("run_param_unc", "Run Parameter Uncertainty Sensitivity", icon = icon("chart-line")),
           plotOutput("param_unc_plot"),
           h4("PRCC (Partial Rank Correlation Coefficient)"),
           actionButton("run_prcc", "Run PRCC Analysis", icon = icon("calculator")),
           plotOutput("prcc_plot"),
           h4("Morris Elementary Effects"),
           fluidRow(
             column(4, numericInput("morris_levels", "Levels", value = 8, min = 4)),
             column(4, numericInput("morris_traj", "Trajectories", value = 30, min = 10)),
             column(4, checkboxInput("morris_global", "Include dose‑response parameters", value = FALSE))
           ),
           actionButton("run_morris", "Run Morris Analysis", icon = icon("diagram-project")),
           plotOutput("morris_plot"),
           hr(),
           h4("Global Sensitivity (Variability + Parameter)"),
           p("Requires 2D simulation with recorded inputs. Click 'Prepare Data' first, then run analyses."),
           actionButton("run_global_prep", "Prepare Global Sensitivity Data", icon = icon("database")),
           br(), br(),
           actionButton("run_global_spearman", "Global Spearman", icon = icon("chart-bar")),
           plotOutput("global_spearman_plot"),
           br(),
           actionButton("run_global_prcc", "Global PRCC", icon = icon("chart-bar")),
           plotOutput("global_prcc_plot"),
           br(),
           actionButton("run_global_morris", "Global Morris", icon = icon("diagram-project")),
           plotOutput("global_morris_plot"),
           h4("Sobol' Variance Decomposition (requires 2D simulation)"),
           fluidRow(
             column(4, numericInput("sobol_n", "Sobol' Sample Size (n)", value = 1000, min = 500, step = 500)),
             column(4, numericInput("sobol_nboot", "Bootstrap Replicates", value = 100, min = 50, step = 50))
           ),
           plotOutput("sobol_plot"),
           verbatimTextOutput("sobol_text")
  ),

  tabPanel("Scenario Comparison",
           sidebarLayout(
             sidebarPanel(
               h4("Scenario Matrix"),
               checkboxGroupInput("comp_groups", "Population Groups:",
                                  choices = c("Adult", "Male", "Female", "Child", "18Y-44Y", "Urban", "Rural"),
                                  selected = c("Adult")),
               checkboxGroupInput("comp_micro", "Microenvironments:",
                                  choices = c("Indoor", "Outdoor"),
                                  selected = "Indoor"),
               numericInput("conc_multiplier", "Concentration Multiplier (x current)",
                            value = 1, min = 0.1, max = 100, step = 0.1),
               hr(),
               numericInput("comp_n_iter", "1D Iterations per Scenario", value = 2000, min = 500),
               actionButton("run_comparison", "Run Comparison", icon = icon("chart-bar"), class = "btn-primary")
             ),
             mainPanel(
               plotOutput("comparison_plot"),
               verbatimTextOutput("comparison_summary")
             )
           )
  )
)

# Server
server <- function(input, output, session) {

  # ---------- Custom publication-ready tornado plots ----------
  custom_tornado <- function(ranking, title, xlab = "Correlation") {
    val_col <- if ("PRCC" %in% names(ranking)) "PRCC" else "Correlation"
    ranking <- ranking[!is.na(ranking[[val_col]]), , drop = FALSE]
    if (nrow(ranking) == 0) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "All correlations are NA") + theme_void())
    }
    ggplot(ranking, aes(x = reorder(Variable, .data[[val_col]]),
                        y = .data[[val_col]])) +
      geom_col(fill = ifelse(ranking[[val_col]] > 0, "#2c7bb6", "#d7191c"),
               width = 0.6) +
      coord_flip() +
      scale_y_continuous(expand = expansion(mult = c(0.006, 0.05))) +
      labs(title = title, x = "", y = xlab) +
      theme_classic(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        axis.text = element_text(color = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()
      )
  }

  custom_tornado_log <- function(ee_df, title, xlab = expression(mu * "*")) {
    df <- data.frame(Variable = rownames(ee_df), mu_star = ee_df[, "mu.star"])
    df <- df[!is.na(df$mu_star), , drop = FALSE]
    if (nrow(df) == 0) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "All mu* are NA") + theme_void())
    }
    ggplot(df, aes(x = reorder(Variable, mu_star), y = mu_star)) +
      geom_col(fill = "#5e3c99", width = 0.6) +
      scale_y_log10(
        expand = expansion(mult = c(0.006, 0.05)),
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x)),
        minor_breaks = NULL
      ) +
      annotation_logticks(base = 10, sides = "b", long = unit(0.15, "cm"),
                          mid = unit(0.1, "cm"), short = unit(0.05, "cm")) +
      coord_flip() +
      labs(title = title, x = "", y = xlab) +
      theme_classic(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        axis.text = element_text(color = "black"),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()
      )
  }

  # Dynamically load population groups
  observe({
    all_groups <- sort(unique(airqmra::ir_long_table$classification))
    updateSelectInput(session, "exposure_group", choices = all_groups, selected = "Adult")
  })

  # Dynamically update column names of uploaded file
  observeEvent(input$conc_file, {
    header <- read.csv(input$conc_file$datapath, nrows = 0)
    updateSelectInput(session, "conc_col", choices = names(header))
  })

  # Truncated survival distribution (ensures values in [0,1])
  surv_dist_trunc <- reactive({
    surv_rec <- get_survival(input$pathogen)
    if (surv_rec$dist == "lognormal") {
      function(n) {
        pmax(0, pmin(1, rlnorm(n, meanlog = surv_rec$log_mean, sdlog = surv_rec$log_sd)))
      }
    } else if (surv_rec$dist == "constant") {
      function(n) rep(surv_rec$value, n)
    } else {
      function(n) rep(0.5, n)
    }
  })

  # Robust concentration distribution function
  conc_dist <- reactive({
    req(input$conc_type)
    if (input$conc_type == "default") {
      function(n) rlnorm(n, meanlog = log(50), sdlog = log(1.8))
    } else if (input$conc_type == "fixed") {
      req(input$conc_fixed)
      function(n) rep(input$conc_fixed, n)
    } else if (input$conc_type == "file") {
      req(input$conc_file, input$conc_col, input$conc_method)
      fpath <- input$conc_file$datapath
      if (is.null(fpath) || !file.exists(fpath) || file.size(fpath) == 0) {
        showNotification("Uploaded file is empty or missing. Using default concentration.", type = "warning")
        return(function(n) rlnorm(n, meanlog = log(50), sdlog = log(1.8)))
      }
      tryCatch(
        import_concentration(fpath, col = input$conc_col, method = input$conc_method),
        error = function(e) {
          showNotification(paste("Concentration import failed:", e$message), type = "error")
          function(n) rlnorm(n, meanlog = log(50), sdlog = log(1.8))
        }
      )
    }
  })

  # Build scenario (used everywhere except when explicitly overridden)
  build_scenario <- reactive({
    scen <- new_scenario(
      pathogen_id          = input$pathogen,
      exposure_group       = input$exposure_group,
      microenvironment     = input$microenvironment,
      season               = input$season,
      population_deposition = input$population_deposition,
      use_uncertainty      = as.logical(input$use_uncertainty),
      C_air_dist           = conc_dist(),
      survival             = surv_dist_trunc(),
      dr_generator         = NULL
    )
    scen
  })

  # ---- Helper: extract model type and parameter names ----
  dose_response_info <- reactive({
    rec <- pathogen_info(input$pathogen)
    if (nrow(rec) == 0) return(list(model = "unknown", params = ""))
    model <- rec$model[1]
    if (model == "exponential") {
      list(model = model, params = "k")
    } else if (model %in% c("exact_beta_poisson", "beta_poisson")) {
      list(model = model, params = "alpha, beta")
    } else if (model == "approx_beta_poisson") {
      list(model = model, params = "alpha, N50")
    } else {
      list(model = model, params = "")
    }
  })

  output$scenario_summary <- renderPrint({
    scen <- build_scenario()
    cat("Pathogen:", scen$pathogen, "\n")
    model_type <- NULL
    params <- NULL
    if (!is.null(scen$dr_pars)) {
      model_type <- scen$dr_pars$type
      params <- scen$dr_pars
    } else if (!is.null(scen$pathogen_id)) {
      rec <- pathogen_info(scen$pathogen_id)
      if (nrow(rec) > 0) {
        model_type <- rec$model[1]
        if (model_type == "exponential") {
          params <- list(k = rec$k[1])
        } else if (model_type == "exact_beta_poisson") {
          params <- list(alpha = rec$alpha[1], beta = rec$beta[1])
        } else if (model_type == "approx_beta_poisson") {
          params <- list(alpha = rec$alpha[1], N50 = rec$N50[1])
        }
      }
    }
    cat("Model type:", ifelse(is.null(model_type), "unknown", model_type), "\n")
    if (!is.null(params)) {
      if (model_type == "exponential") {
        k_val <- params$k
        if (is.numeric(k_val) && is.finite(k_val)) cat("k =", signif(k_val, 3), "\n")
      } else if (model_type == "exact_beta_poisson") {
        alpha_val <- params$alpha; beta_val <- params$beta
        if (is.numeric(alpha_val) && is.numeric(beta_val))
          cat("alpha =", alpha_val, " beta =", signif(beta_val, 3), "\n")
      } else if (model_type == "approx_beta_poisson") {
        alpha_val <- params$alpha; n50_val <- params$N50
        if (is.numeric(alpha_val) && is.numeric(n50_val)) {
          n50_str <- if (is.finite(n50_val)) signif(n50_val, 3) else "unknown"
          cat("alpha =", alpha_val, " N50 =", n50_str, "\n")
        }
      }
    } else {
      cat("Uncertainty enabled (parameters sampled from distributions).\n")
    }
    if (input$conc_type == "file") {
      cat("\nImported concentration distribution info:\n")
      print(conc_dist())
    }
  })

  output$scenario_check <- renderPrint({
    scen <- build_scenario()
    cat("Running scenario checks...\n")
    check_scenario(scen)
  })

  # 1D simulation
  res1d <- eventReactive(input$run1d, {
    scen <- build_scenario()
    single <- run_qmra_1d(scen, n_iter = input$n_iter)
    if (input$annual_type == "fixed_freq") {
      freq <- input$freq
    } else {
      freq <- function(n) rpois(n, lambda = input$lambda)
    }
    annual <- annual_risk(single, frequency = freq)
    list(single = single, annual = annual, scen = scen)
  })

  output$hist1d <- renderPlot({
    req(res1d())
    ggplot(data.frame(risk = res1d()$single), aes(x = risk)) +
      geom_density(fill = "steelblue", alpha = 0.4, color = "steelblue") +
      scale_x_log10() +
      labs(title = "Single-Exposure Infection Risk Density",
           x = "Infection Probability (log scale)", y = "Density") +
      theme_light(base_size = 14)
  })

  output$summary1d <- renderPrint({
    req(res1d())
    summary(res1d()$single)
  })

  output$annual_hist <- renderPlot({
    req(res1d())
    ggplot(data.frame(risk = res1d()$annual), aes(x = risk)) +
      geom_density(fill = "forestgreen", alpha = 0.4, color = "forestgreen") +
      scale_x_log10() +
      labs(title = "Annual Infection Risk Density",
           x = "Annual Infection Probability (log scale)", y = "Density") +
      theme_light(base_size = 14)
  })

  output$annual_summary <- renderPrint({
    req(res1d())
    summary(res1d()$annual)
  })

  # 2D simulation
  res2d <- eventReactive(input$run2d, {
    scen_unc <- new_scenario(
      pathogen_id          = input$pathogen,
      exposure_group       = input$exposure_group,
      microenvironment     = input$microenvironment,
      season               = input$season,
      population_deposition = input$population_deposition,
      use_uncertainty      = TRUE,
      C_air_dist           = conc_dist(),
      survival             = surv_dist_trunc(),
      dr_generator         = NULL
    )
    if (input$annual_type == "fixed_freq") {
      freq <- input$freq
    } else {
      freq <- function(n) rpois(n, lambda = input$lambda)
    }
    withProgress(message = "Running 2D simulation...", value = 0, {
      res <- run_qmra_2d(scen_unc, n_var = input$n_var, n_unc = input$n_unc, frequency = freq)
    })
    res
  })

  output$hist2d_overall <- renderPlot({
    req(res2d())
    ggplot(data.frame(risk = res2d()$overall_risk), aes(x = risk)) +
      geom_density(fill = "orange", alpha = 0.4, color = "orange") +
      scale_x_log10() +
      labs(title = "Combined Risk Distribution (Variability + Uncertainty)",
           x = "Infection Probability (log scale)", y = "Density") +
      theme_light(base_size = 14)
  })

  output$summary2d <- renderPrint({
    req(res2d())
    summary(res2d()$overall_risk)
  })

  output$hist2d_mean <- renderPlot({
    req(res2d())
    ggplot(data.frame(mean_risk = res2d()$mean_risk), aes(x = mean_risk)) +
      geom_density(fill = "purple", alpha = 0.4, color = "purple") +
      scale_x_log10() +
      labs(title = "Uncertainty Distribution of Mean Risk",
           x = "Mean Infection Probability (log scale)", y = "Density") +
      theme_light(base_size = 14)
  })

  output$var_boxplot <- renderPlot({
    req(res2d())
    mat <- res2d()$risk_matrix
    ncols <- ncol(mat)
    sample_cols <- if (ncols <= 20) 1:ncols else sort(sample.int(ncols, 20))
    df <- as.data.frame(mat[, sample_cols, drop = FALSE])
    colnames(df) <- paste0("Iter ", sample_cols)
    df_long <- tidyr::pivot_longer(df, everything(), names_to = "Iteration", values_to = "Risk")

    ggplot(df_long, aes(x = Iteration, y = Risk)) +
      geom_boxplot(outlier.size = 0.5, fill = "#87CEEB", alpha = 0.7) +
      scale_y_log10() +
      labs(title = "Variability of Risk across Uncertainty Iterations",
           subtitle = paste0("Randomly selected ", ncol(df), " iterations"),
           x = "Uncertainty Iteration", y = "Infection Probability (log scale)") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })

  # ---- Global Sensitivity Data Preparation ----
  global_res2d <- eventReactive(input$run_global_prep, {
    scen_unc <- new_scenario(
      pathogen_id          = input$pathogen,
      exposure_group       = input$exposure_group,
      microenvironment     = input$microenvironment,
      season               = input$season,
      population_deposition = input$population_deposition,
      use_uncertainty      = TRUE,
      C_air_dist           = conc_dist(),
      survival             = surv_dist_trunc(),
      dr_generator         = NULL
    )
    freq <- if (input$annual_type == "fixed_freq") input$freq else 20
    showNotification("Running small 2D simulation to record inputs...", type = "default")
    withProgress(message = "Preparing global sensitivity data", value = 0, {
      res <- run_qmra_2d(scen_unc, n_var = 500, n_unc = 200, frequency = freq,
                         record_inputs = TRUE)
    })
    res
  })

  # ---- Global Spearman ----
  observeEvent(input$run_global_spearman, {
    req(global_res2d())
    sens <- sensitivity_global(global_res2d(), method = "spearman")
    dr_info <- dose_response_info()
    output$global_spearman_plot <- renderPlot({
      custom_tornado(sens$ranking,
                     title = paste0("Global Spearman (incl. ", dr_info$params, ")"),
                     xlab = "Spearman correlation")
    })
  })

  # ---- Global PRCC ----
  observeEvent(input$run_global_prcc, {
    req(global_res2d())
    sens <- sensitivity_global_prcc(global_res2d(), nboot = 200)
    dr_info <- dose_response_info()
    output$global_prcc_plot <- renderPlot({
      custom_tornado(sens$ranking,
                     title = paste0("Global PRCC (incl. ", dr_info$params, ")"),
                     xlab = "Partial rank correlation")
    })
  })

  # ---- Global Morris ----
  observeEvent(input$run_global_morris, {
    scen_unc <- new_scenario(
      pathogen_id          = input$pathogen,
      exposure_group       = input$exposure_group,
      microenvironment     = input$microenvironment,
      season               = input$season,
      population_deposition = input$population_deposition,
      use_uncertainty      = TRUE,
      C_air_dist           = conc_dist(),
      survival             = surv_dist_trunc(),
      dr_generator         = NULL
    )
    sens <- sensitivity_global_morris(scen_unc, n_levels = 8, n_trajectories = 50)
    dr_info <- dose_response_info()
    output$global_morris_plot <- renderPlot({
      custom_tornado_log(sens$ee,
                         title = paste0("Global Morris (mu*, incl. ", dr_info$params, ")"))
    })
  })

  # Sensitivity analysis - Spearman
  observe({
    output$tornado1d <- renderPlot({
      scen_point <- new_scenario(
        pathogen_id          = input$pathogen,
        exposure_group       = input$exposure_group,
        microenvironment     = input$microenvironment,
        season               = input$season,
        population_deposition = input$population_deposition,
        use_uncertainty      = FALSE,
        C_air_dist           = conc_dist(),
        survival             = surv_dist_trunc(),
        dr_generator         = NULL
      )
      sens <- sensitivity_1d(scen_point, n_iter = 2000)
      custom_tornado(sens$sensitivity$ranking,
                     title = "1D Spearman Sensitivity",
                     xlab = "Spearman correlation")
    })
  })

  # Parameter Uncertainty – multi‑parameter version
  param_unc_all <- eventReactive(input$run_param_unc, {
    req(res2d())
    param_list <- res2d()$param_list
    mean_risk <- res2d()$mean_risk
    # Gather all unique parameter names
    param_names <- unique(unlist(lapply(param_list, names)))
    cor_list <- lapply(param_names, function(pn) {
      vals <- sapply(param_list, function(x) x[pn])
      cor(vals, mean_risk, method = "spearman", use = "complete.obs")
    })
    cors <- unlist(cor_list)
    ranking <- data.frame(
      Variable = param_names,
      Correlation = cors,
      AbsCorrelation = abs(cors),
      stringsAsFactors = FALSE
    )
    ranking <- ranking[order(-ranking$AbsCorrelation), ]
    ranking
  })

  output$param_unc_plot <- renderPlot({
    req(param_unc_all())
    custom_tornado(param_unc_all(),
                   title = "Parameter Uncertainty (all parameters)",
                   xlab = "Spearman correlation")
  })

  # PRCC analysis
  prcc_result <- eventReactive(input$run_prcc, {
    scen_point <- new_scenario(
      pathogen_id          = input$pathogen,
      exposure_group       = input$exposure_group,
      microenvironment     = input$microenvironment,
      season               = input$season,
      population_deposition = input$population_deposition,
      use_uncertainty      = FALSE,
      C_air_dist           = conc_dist(),
      survival             = surv_dist_trunc(),
      dr_generator         = NULL
    )
    sensitivity_prcc(scen_point, n_iter = 3000, nboot = 200)
  })

  output$prcc_plot <- renderPlot({
    req(prcc_result())
    custom_tornado(prcc_result()$ranking,
                   title = "PRCC Sensitivity",
                   xlab = "Partial rank correlation")
  })

  # Morris analysis (global or local)
  morris_result <- eventReactive(input$run_morris, {
    if (input$morris_global) {
      scen_unc <- new_scenario(
        pathogen_id          = input$pathogen,
        exposure_group       = input$exposure_group,
        microenvironment     = input$microenvironment,
        season               = input$season,
        population_deposition = input$population_deposition,
        use_uncertainty      = TRUE,
        C_air_dist           = conc_dist(),
        survival             = surv_dist_trunc(),
        dr_generator         = NULL
      )
      sensitivity_global_morris(scen_unc, n_levels = input$morris_levels,
                                n_trajectories = input$morris_traj)
    } else {
      scen_point <- new_scenario(
        pathogen_id          = input$pathogen,
        exposure_group       = input$exposure_group,
        microenvironment     = input$microenvironment,
        season               = input$season,
        population_deposition = input$population_deposition,
        use_uncertainty      = FALSE,
        C_air_dist           = conc_dist(),
        survival             = surv_dist_trunc(),
        dr_generator         = NULL
      )
      sensitivity_morris(scen_point, n_levels = input$morris_levels,
                         n_trajectories = input$morris_traj)
    }
  })

  output$morris_plot <- renderPlot({
    req(morris_result())
    custom_tornado_log(morris_result()$ee,
                       title = if (input$morris_global) "Morris Sensitivity (mu*, incl. parameter)" else "Morris Sensitivity (mu*)")
  })

  # Sobol' analysis
  sobol_result <- reactive({
    req(res2d(), input$sobol_n, input$sobol_nboot)
    tryCatch({
      scen <- new_scenario(
        pathogen_id          = input$pathogen,
        exposure_group       = input$exposure_group,
        microenvironment     = input$microenvironment,
        season               = input$season,
        population_deposition = input$population_deposition,
        use_uncertainty      = TRUE,
        C_air_dist           = conc_dist(),
        dr_generator         = NULL,
        survival             = surv_dist_trunc()
      )
      sobol_sensitivity(scen, n = input$sobol_n, nboot = input$sobol_nboot)
    }, error = function(e) {
      showNotification(paste("Sobol' analysis failed:", e$message), type = "error")
      NULL
    })
  })

  output$sobol_plot <- renderPlot({
    req(sobol_result())
    res <- sobol_result()
    S <- as.data.frame(res$S) %>% rownames_to_column("Variable") %>%
      rename(original = original)
    T <- as.data.frame(res$T) %>% rownames_to_column("Variable") %>%
      rename(original = original)
    has_ci <- all(c("min. c.i.", "max. c.i.") %in% names(S))
    if (has_ci) {
      S <- S %>% rename(min_ci = `min. c.i.`, max_ci = `max. c.i.`)
      T <- T %>% rename(min_ci = `min. c.i.`, max_ci = `max. c.i.`)
    }
    combined <- bind_rows(S %>% mutate(Type = "First-order"),
                          T %>% mutate(Type = "Total-effect"))
    p <- ggplot(combined, aes(x = reorder(Variable, original), y = original, fill = Type)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6) +
      coord_flip() +
      scale_fill_manual(values = c("First-order" = "steelblue", "Total-effect" = "coral")) +
      labs(title = "Sobol' Variance Decomposition",
           subtitle = if (has_ci) "First-order and Total-effect Indices (95% CI)" else "First-order and Total-effect Indices (no bootstrap CI)",
           x = "", y = "Sobol' Index") +
      theme_minimal(base_size = 14)
    if (has_ci) {
      p <- p + geom_errorbar(aes(ymin = pmax(0, min_ci), ymax = max_ci),
                             position = position_dodge(0.7), width = 0.15, color = "grey30")
    }
    p
  })

  output$sobol_text <- renderPrint({
    req(sobol_result())
    print(sobol_result())
  })

  # Scenario comparison
  comparison_result <- eventReactive(input$run_comparison, {
    base_scen <- build_scenario()
    pathogen_id <- base_scen$pathogen_id
    season      <- base_scen$season
    pop_dep     <- base_scen$population_deposition

    orig_conc_fun <- conc_dist()
    multiplier    <- input$conc_multiplier
    conc_fun <- if (multiplier == 1) orig_conc_fun else {
      function(n) orig_conc_fun(n) * multiplier
    }

    groups <- input$comp_groups
    micros <- input$comp_micro
    if (length(groups) == 0 || length(micros) == 0) return(NULL)

    grid <- expand.grid(group = groups, micro = micros, stringsAsFactors = FALSE)
    scens <- list()
    for (i in seq_len(nrow(grid))) {
      scens[[paste0(grid$group[i], ", ", grid$micro[i])]] <- new_scenario(
        pathogen_id          = pathogen_id,
        exposure_group       = grid$group[i],
        microenvironment     = grid$micro[i],
        season               = season,
        population_deposition = pop_dep,
        use_uncertainty      = FALSE,
        C_air_dist           = conc_fun,
        survival             = surv_dist_trunc(),
        dr_generator         = NULL
      )
    }
    if (length(scens) == 0) return(NULL)
    batch_simulate_1d(scens, n_iter = input$comp_n_iter, frequency = 20)
  })

  output$comparison_plot <- renderPlot({
    req(comparison_result())
    plot_scenario_comparison(comparison_result()$raw, plot_type = "box", log_scale = TRUE) +
      ggtitle("Single Exposure Risk Comparison")
  })

  output$comparison_summary <- renderPrint({
    req(comparison_result())
    comparison_result()$summary
  })

  # Set output directory
  observeEvent(input$set_dir, {
    updateTextInput(session, "report_dir", value = getwd())
  })

  # Report generation
  observeEvent(input$generate_report, {
    if (is.null(res1d())) {
      showNotification("Please run the 1D Simulation first.", type = "error")
      return()
    }

    scen_pt <- new_scenario(
      pathogen_id          = input$pathogen,
      exposure_group       = input$exposure_group,
      microenvironment     = input$microenvironment,
      season               = input$season,
      population_deposition = input$population_deposition,
      use_uncertainty      = FALSE,
      C_air_dist           = conc_dist(),
      survival             = surv_dist_trunc(),
      dr_generator         = NULL
    )

    results <- tryCatch({
      core <- list(
        single_risk    = res1d()$single,
        annual_risk    = res1d()$annual,
        sensitivity_1d = sensitivity_1d(scen_pt, n_iter = 2000)
      )
      if (!is.null(prcc_result())) core$prcc <- prcc_result()
      if (!is.null(morris_result())) core$morris <- morris_result()
      if (!is.null(sobol_result())) core$sobol <- sobol_result()
      if (!is.null(res2d())) core$res2d <- res2d()
      core
    }, error = function(e) {
      showNotification(paste("Failed to prepare report data:", e$message), type = "error")
      return(NULL)
    })
    if (is.null(results)) return()

    out_dir <- input$report_dir
    if (!dir.exists(out_dir)) {
      out_dir <- getwd()
    }
    output_file <- file.path(out_dir, paste0("airqmra_report_", Sys.Date()))
    if (input$report_format == "word") {
      output_file <- paste0(output_file, ".docx")
    } else {
      output_file <- paste0(output_file, ".html")
    }

    showNotification("Generating report...", type = "default")
    withProgress(message = "Rendering report", value = 0, {
      success <- FALSE
      tryCatch({
        generate_risk_report(
          scenario = scen_pt,
          risk_results = results,
          output_file = output_file,
          output_format = input$report_format,
          language = input$report_lang
        )
        success <- TRUE
      }, error = function(e) {
        showNotification(paste("Report generation failed:", e$message), type = "error", duration = 10)
      })
    })

    if (success) {
      showNotification(paste("Report saved to:", output_file), type = "default", duration = 10)
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)
