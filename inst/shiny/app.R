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
               p("2D simulation is time-consuming. Please be patient.")
             ),
             mainPanel(
               h4("Current Scenario Summary"),
               verbatimTextOutput("scenario_summary")
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
           h4("Sobol' Variance Decomposition (requires 2D simulation)"),
           fluidRow(
             column(4, numericInput("sobol_n", "Sobol' Sample Size (n)", value = 1000, min = 500, step = 500)),
             column(4, numericInput("sobol_nboot", "Bootstrap Replicates", value = 100, min = 50, step = 50))
           ),
           plotOutput("sobol_plot"),
           verbatimTextOutput("sobol_text")
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

  # 2D simulation (uses user-specified n_var and n_unc)
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

  # Sensitivity analysis
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

  # Sobol' analysis
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
}

# Run the application
shinyApp(ui = ui, server = server)
