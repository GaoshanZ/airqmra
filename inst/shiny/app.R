library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(airqmra)

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
               radioButtons("use_uncertainty", "Dose‑Response Parameter Uncertainty",
                            choices = c("Point Estimate" = FALSE, "Include Uncertainty" = TRUE),
                            selected = FALSE),
               hr(),
               h5("Air Concentration (copies/m³)"),
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
           h4("PRCC (Partial Rank Correlation Coefficient)"),
           actionButton("run_prcc", "Run PRCC Analysis", icon = icon("calculator")),
           plotOutput("prcc_plot"),
           h4("Morris Elementary Effects"),
           fluidRow(
             column(4, numericInput("morris_levels", "Levels", value = 8, min = 4)),
             column(4, numericInput("morris_traj", "Trajectories", value = 30, min = 10))
           ),
           actionButton("run_morris", "Run Morris Analysis", icon = icon("diagram-project")),
           plotOutput("morris_plot"),
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
               numericInput("conc_multiplier", "Concentration Multiplier (× current)",
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

  # Dynamically load population groups
  observe({
    all_groups <- sort(unique(airqmra::ir_long_table$classification))
    updateSelectInput(session, "exposure_group", choices = all_groups)
  })

  # Dynamically update column names of uploaded file
  observeEvent(input$conc_file, {
    header <- read.csv(input$conc_file$datapath, nrows = 0)
    updateSelectInput(session, "conc_col", choices = names(header))
  })

  # Concentration distribution function
  conc_dist <- reactive({
    req(input$conc_type)
    if (input$conc_type == "default") {
      function(n) rlnorm(n, meanlog = log(50), sdlog = log(1.8))
    } else if (input$conc_type == "fixed") {
      req(input$conc_fixed)
      function(n) rep(input$conc_fixed, n)
    } else if (input$conc_type == "file") {
      req(input$conc_file, input$conc_col, input$conc_method)
      import_concentration(input$conc_file$datapath,
                           col = input$conc_col,
                           method = input$conc_method)
    }
  })

  # Build scenario
  build_scenario <- reactive({
    scen <- new_scenario(
      pathogen_id          = input$pathogen,
      exposure_group       = input$exposure_group,
      microenvironment     = input$microenvironment,
      season               = input$season,
      population_deposition = input$population_deposition,
      use_uncertainty      = as.logical(input$use_uncertainty),
      C_air_dist           = conc_dist()
    )
    scen
  })

  output$scenario_summary <- renderPrint({
    scen <- build_scenario()
    cat("Pathogen:", scen$pathogen, "\n")
    cat("Model type:", scen$dr_pars$type, "\n")
    if (scen$dr_pars$type == "exponential") {
      cat("k =", signif(scen$dr_pars$k, 3), "\n")
    } else if (scen$dr_pars$type == "exact_beta_poisson") {
      cat("alpha =", scen$dr_pars$alpha, " beta =", signif(scen$dr_pars$beta, 3), "\n")
    } else if (scen$dr_pars$type == "approx_beta_poisson") {
      cat("alpha =", scen$dr_pars$alpha, " N50 =", signif(scen$dr_pars$N50, 3), "\n")
    }
    if (input$conc_type == "file") {
      cat("\nImported concentration distribution info:\n")
      print(conc_dist())
    }
  })

  # Scenario validation
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
      labs(title = "Single‑Exposure Infection Risk Density",
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
      C_air_dist           = conc_dist()
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

  # Variability boxplot across uncertainty iterations
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

  # Sensitivity analysis - Spearman (always available)
  observe({
    output$tornado1d <- renderPlot({
      scen_point <- new_scenario(
        pathogen_id          = input$pathogen,
        exposure_group       = input$exposure_group,
        microenvironment     = input$microenvironment,
        season               = input$season,
        population_deposition = input$population_deposition,
        use_uncertainty      = FALSE,
        C_air_dist           = conc_dist()
      )
      sens <- sensitivity_1d(scen_point, n_iter = 2000)
      tornado_plot(sens, title = "1D Spearman Sensitivity")
    })
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
      C_air_dist           = conc_dist()
    )
    sensitivity_prcc(scen_point, n_iter = 3000, nboot = 200)
  })

  output$prcc_plot <- renderPlot({
    req(prcc_result())
    tornado_plot(prcc_result(), title = "PRCC Sensitivity")
  })

  # Morris analysis
  morris_result <- eventReactive(input$run_morris, {
    scen_point <- new_scenario(
      pathogen_id          = input$pathogen,
      exposure_group       = input$exposure_group,
      microenvironment     = input$microenvironment,
      season               = input$season,
      population_deposition = input$population_deposition,
      use_uncertainty      = FALSE,
      C_air_dist           = conc_dist()
    )
    sensitivity_morris(scen_point, n_levels = input$morris_levels,
                       n_trajectories = input$morris_traj)
  })

  output$morris_plot <- renderPlot({
    req(morris_result())
    tornado_plot(morris_result(), title = "Morris Sensitivity (mu*)", log_scale = TRUE)
  })

  # Sobol' analysis (unchanged)
  sobol_result <- reactive({
    req(res2d(), input$sobol_n, input$sobol_nboot)
    scen <- new_scenario(
      pathogen_id          = input$pathogen,
      exposure_group       = input$exposure_group,
      microenvironment     = input$microenvironment,
      season               = input$season,
      population_deposition = input$population_deposition,
      use_uncertainty      = TRUE,
      C_air_dist           = conc_dist()
    )
    sobol_sensitivity(scen, n = input$sobol_n, nboot = input$sobol_nboot)
  })

  output$sobol_plot <- renderPlot({
    req(sobol_result())
    res <- sobol_result()
    S <- as.data.frame(res$S) %>% rownames_to_column("Variable") %>%
      rename(original = original, min_ci = `min. c.i.`, max_ci = `max. c.i.`)
    T <- as.data.frame(res$T) %>% rownames_to_column("Variable") %>%
      rename(original = original, min_ci = `min. c.i.`, max_ci = `max. c.i.`)
    combined <- bind_rows(S %>% mutate(Type = "First-order"),
                          T %>% mutate(Type = "Total-effect"))
    ggplot(combined, aes(x = reorder(Variable, original), y = original, fill = Type)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6) +
      geom_errorbar(aes(ymin = pmax(0, min_ci), ymax = max_ci),
                    position = position_dodge(width = 0.7), width = 0.15, color = "grey30") +
      coord_flip() +
      scale_fill_manual(values = c("First-order" = "steelblue", "Total-effect" = "coral")) +
      labs(title = "Sobol' Variance Decomposition",
           subtitle = "First-order and Total-effect Indices (95% CI)",
           x = "", y = "Sobol' Index") +
      theme_minimal(base_size = 14)
  })

  output$sobol_text <- renderPrint({
    req(sobol_result())
    print(sobol_result())
  })

  # Scenario comparison (batch 1D)
  # Scenario comparison: adaptive matrix based on user selections
  comparison_result <- eventReactive(input$run_comparison, {
    # Retrieve fixed parameters from the current scenario (pathogen, season, deposition)
    base_scen <- build_scenario()
    pathogen_id <- base_scen$pathogen_id
    season      <- base_scen$season
    pop_dep     <- base_scen$population_deposition

    # Use the current concentration distribution (with optional multiplier)
    orig_conc_fun <- conc_dist()
    multiplier    <- input$conc_multiplier
    conc_fun <- if (multiplier == 1) orig_conc_fun else {
      function(n) orig_conc_fun(n) * multiplier
    }

    groups <- input$comp_groups
    micros <- input$comp_micro
    if (length(groups) == 0 || length(micros) == 0) return(NULL)

    # Generate all combinations of selected groups and microenvironments
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
        C_air_dist           = conc_fun
      )
    }
    if (length(scens) == 0) return(NULL)

    # Run batch 1D simulation with user-defined iterations and annual frequency
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

  # Report generation
  observeEvent(input$generate_report, {
    scen <- build_scenario()
    # Collect available results
    results <- list(single_risk = NULL, annual_risk = NULL,
                    sensitivity_1d = NULL, prcc = NULL, morris = NULL,
                    sobol = NULL, res2d = NULL)
    if (!is.null(res1d())) {
      results$single_risk <- res1d()$single
      results$annual_risk <- res1d()$annual
    }
    if (!is.null(tornado1d())) { } # tornado1d is reactive not a stored value, we can recalc
    # Simplified: we recalc sens1d if needed, but to save time, we just pass what we have
    # For report generation, we can run a small 1D if not present
    if (is.null(results$single_risk)) {
      showNotification("Running a quick 1D simulation for the report...", type = "message")
      s <- run_qmra_1d(scen, n_iter = 2000)
      results$single_risk <- s
      results$annual_risk <- annual_risk(s, frequency = 20)
    }
    if (is.null(results$sensitivity_1d)) {
      results$sensitivity_1d <- sensitivity_1d(scen, n_iter = 2000)
    }
    # Optional: add prcc, morris, sobol if already computed
    if (!is.null(prcc_result())) results$prcc <- prcc_result()
    if (!is.null(morris_result())) results$morris <- morris_result()
    if (!is.null(sobol_result())) results$sobol <- sobol_result()
    if (!is.null(res2d())) results$res2d <- res2d()

    output_file <- paste0("airqmra_report_", Sys.Date(), ".", input$report_format)
    if (input$report_format == "word") output_file <- paste0("airqmra_report_", Sys.Date(), ".docx")

    generate_risk_report(
      scenario = scen,
      risk_results = results,
      output_file = output_file,
      output_format = input$report_format,
      language = input$report_lang
    )
    showNotification(paste("Report generated:", output_file), type = "default")
  })

}

# Run the application
shinyApp(ui = ui, server = server)
