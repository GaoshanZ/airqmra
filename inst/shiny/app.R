library(shiny)
library(ggplot2)
library(airqmra)

# 病原体列表（从数据库中提取）
pathogen_choices <- unique(airqmra::pathogen_db$pathogen_id)
names(pathogen_choices) <- paste0(airqmra::pathogen_db$pathogen_cn[match(pathogen_choices, airqmra::pathogen_db$pathogen_id)],
                                  " (", pathogen_choices, ")")

# UI
ui <- navbarPage(
  title = "airqmra 吸入风险评估",
  theme = shinythemes::shinytheme("flatly"),

  # 选项卡1：场景设置
  tabPanel("场景设置",
           sidebarLayout(
             sidebarPanel(
               selectInput("pathogen", "病原体", choices = pathogen_choices),
               selectInput("exposure_group", "人群分组",
                           choices = c("Adult", "Male", "Female", "18Y-44Y", "Urban", "Rural")),
               selectInput("microenvironment", "微环境",
                           choices = c("Indoor", "Outdoor")),
               selectInput("season", "季节",
                           choices = c("Annual", "Spring & Autumn", "Summer", "Winter")),
               selectInput("population_deposition", "沉积参数人群",
                           choices = c("Male", "Female", "Child", "Senior")),
               radioButtons("use_uncertainty", "剂量-反应参数不确定性",
                            choices = c("点估计" = FALSE, "纳入不确定性" = TRUE),
                            selected = FALSE),
               hr(),
               h5("暴露浓度 (copies/m³)"),
               radioButtons("conc_type", "浓度类型",
                            choices = c("默认对数正态" = "default",
                                        "固定值" = "fixed",
                                        "上传文件" = "file"),
                            selected = "default"),
               conditionalPanel(
                 condition = "input.conc_type == 'fixed'",
                 numericInput("conc_fixed", "固定浓度值", value = 50, min = 0)
               ),
               conditionalPanel(
                 condition = "input.conc_type == 'file'",
                 fileInput("conc_file", "选择 CSV 文件",
                           accept = c("text/csv", ".csv")),
                 selectInput("conc_col", "浓度列名", choices = NULL),
                 radioButtons("conc_method", "拟合方法",
                              choices = c("经验分布" = "empirical",
                                          "对数正态" = "lognormal",
                                          "自动选择" = "auto"),
                              selected = "lognormal")
               ),
               hr(),
               h5("年暴露参数"),
               radioButtons("annual_type", "年频率",
                            choices = c("固定次数" = "fixed_freq", "泊松分布" = "poisson"),
                            selected = "fixed_freq"),
               conditionalPanel(
                 condition = "input.annual_type == 'fixed_freq'",
                 numericInput("freq", "每年暴露次数", value = 20, min = 1)
               ),
               conditionalPanel(
                 condition = "input.annual_type == 'poisson'",
                 numericInput("lambda", "泊松均值", value = 20, min = 0)
               ),
               hr(),
               actionButton("run1d", "运行一维模拟", icon = icon("play"), class = "btn-primary"),
               br(), br(),
               actionButton("run2d", "运行二维模拟", icon = icon("random"), class = "btn-success"),
               br(), br(),
               p("二维模拟较耗时，请耐心等待。")
             ),
             mainPanel(
               h4("当前场景摘要"),
               verbatimTextOutput("scenario_summary")
             )
           )
  ),

  # 选项卡2：一维模拟结果
  tabPanel("一维模拟",
           plotOutput("hist1d"),
           verbatimTextOutput("summary1d"),
           h4("年风险"),
           plotOutput("annual_hist"),
           verbatimTextOutput("annual_summary")
  ),

  # 选项卡3：二维模拟结果
  tabPanel("二维模拟",
           plotOutput("hist2d_overall"),
           verbatimTextOutput("summary2d"),
           h4("平均风险的不确定性分布"),
           plotOutput("hist2d_mean")
  ),

  # 选项卡4：敏感性分析
  tabPanel("敏感性分析",
           h4("Spearman 飓风图（一维）"),
           plotOutput("tornado1d"),
           h4("Sobol' 方差分解（需先运行二维模拟）"),
           plotOutput("sobol_plot"),
           verbatimTextOutput("sobol_text")
  )
)

# Server
server <- function(input, output, session) {

  # 动态更新上传文件的列名选择
  observeEvent(input$conc_file, {
    header <- read.csv(input$conc_file$datapath, nrows = 0)
    updateSelectInput(session, "conc_col", choices = names(header))
  })

  # 反应式：浓度分布函数
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

  # 反应式：构建场景对象
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

  # 当前场景摘要
  output$scenario_summary <- renderPrint({
    scen <- build_scenario()
    cat("病原体:", scen$pathogen, "\n")
    cat("模型类型:", scen$dr_pars$type, "\n")
    if (scen$dr_pars$type == "exponential") {
      cat("k =", signif(scen$dr_pars$k, 3), "\n")
    } else if (scen$dr_pars$type == "exact_beta_poisson") {
      cat("alpha =", scen$dr_pars$alpha, " beta =", signif(scen$dr_pars$beta, 3), "\n")
    } else if (scen$dr_pars$type == "approx_beta_poisson") {
      cat("alpha =", scen$dr_pars$alpha, " N50 =", signif(scen$dr_pars$N50, 3), "\n")
    }
    # 如果浓度来自文件，显示拟合信息
    if (input$conc_type == "file") {
      cat("\n导入浓度分布信息:\n")
      print(conc_dist())
    }
  })

  # 一维模拟结果（反应式）
  res1d <- eventReactive(input$run1d, {
    scen <- build_scenario()
    # 确保不确定性参数正确（场景中已经设置）
    single <- run_qmra_1d(scen, n_iter = 5000)

    # 确定频率
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
    hist(res1d()$single, breaks = 50, main = "单次暴露感染风险分布",
         xlab = "感染概率", col = "steelblue")
  })

  output$summary1d <- renderPrint({
    req(res1d())
    summary(res1d()$single)
  })

  output$annual_hist <- renderPlot({
    req(res1d())
    hist(res1d()$annual, breaks = 50, main = "年感染风险分布",
         xlab = "年感染概率", col = "forestgreen")
  })

  output$annual_summary <- renderPrint({
    req(res1d())
    summary(res1d()$annual)
  })

  # 二维模拟结果
  res2d <- eventReactive(input$run2d, {
    # 构建不确定性场景（使用相同的浓度分布）
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

    # 显示进度
    withProgress(message = "正在运行二维模拟...", value = 0, {
      res <- run_qmra_2d(scen_unc, n_var = 1000, n_unc = 300, frequency = freq)
    })
    res
  })

  output$hist2d_overall <- renderPlot({
    req(res2d())
    hist(res2d()$overall_risk, breaks = 50, main = "综合风险分布 (变异性+不确定性)",
         xlab = "感染概率", col = "orange")
  })

  output$summary2d <- renderPrint({
    req(res2d())
    summary(res2d()$overall_risk)
  })

  output$hist2d_mean <- renderPlot({
    req(res2d())
    hist(res2d()$mean_risk, breaks = 30, main = "平均风险的不确定性分布",
         xlab = "平均感染概率", col = "purple")
  })

  # 敏感性分析
  observe({
    # 一维敏感性（点估计模式，使用相同浓度分布）
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
      tornado_plot(sens, title = "一维 Spearman 敏感性")
    })
  })

  # Sobol' 分析（需要二维模拟的不确定性信息）
  sobol_result <- reactive({
    req(res2d())  # 确保二维模拟已运行
    scen <- new_scenario(
      pathogen_id          = input$pathogen,
      exposure_group       = input$exposure_group,
      microenvironment     = input$microenvironment,
      season               = input$season,
      population_deposition = input$population_deposition,
      use_uncertainty      = TRUE,
      C_air_dist           = conc_dist()
    )
    # 使用二维模拟记录的信息或重新运行小的 Sobol'（基于数据库中的不确定性参数）
    sobol_sensitivity(scen, n = 300, order = 2, nboot = 50)
  })

  output$sobol_plot <- renderPlot({
    req(sobol_result())
    plot(sobol_result())
  })

  output$sobol_text <- renderPrint({
    req(sobol_result())
    print(sobol_result())
  })
}

# Run the application
shinyApp(ui = ui, server = server)
