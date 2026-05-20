# 从原始暴露参数 CSV 生成 ir_long_table, ir_short_table, time_table
library(readr)
library(dplyr)
library(tidyr)
library(usethis)

# 读入原始数据
raw <- read_csv("data-raw/exposure_raw.csv",
                col_types = cols(.default = col_character()))

# 数值列转换
numeric_cols <- c("IR_log_mean", "IR_log_sd",
                  "Time_outdoor_log_mean", "Time_outdoor_log_sd",
                  "Time_indoor_log_mean", "Time_indoor_log_sd")
raw <- raw %>% mutate(across(any_of(numeric_cols), as.numeric))

# ================= 1. 长期呼吸速率表 =================
ir_long_table <- raw %>%
  filter(activity == "long") %>%
  select(
    classification = Classification,
    IR_dist = IR_Distribution,
    IR_log_mean,
    IR_log_sd,
    IR_unit,
    ref
  ) %>%
  distinct()

# ================= 2. 短期呼吸速率表 =================
ir_short_table <- raw %>%
  filter(grepl("^Short_", activity)) %>%
  mutate(activity_intensity = sub("^Short_", "", activity)) %>%
  select(
    classification = Classification,
    activity_intensity,
    IR_dist = IR_Distribution,
    IR_log_mean,
    IR_log_sd,
    IR_unit,
    ref
  ) %>%
  distinct()

# ================= 3. 时间活动模式表（修正版）=================
time_long <- raw %>%
  filter(activity == "long") %>%
  select(
    classification = Classification,
    Time_outdoor_log_mean,
    Time_outdoor_log_sd,
    Time_indoor_log_mean,
    Time_indoor_log_sd,
    ref
  ) %>%
  distinct()

# 提取季节信息
time_long <- time_long %>%
  mutate(
    season = case_when(
      classification == "Spring and Autumn" ~ "Spring & Autumn",
      classification == "Summer"            ~ "Summer",
      classification == "Winter"            ~ "Winter",
      TRUE                                  ~ "Annual"
    ),
    # 若 classification 本身是季节名，则标准化为 "Adult"
    classification = ifelse(
      classification %in% c("Spring and Autumn", "Summer", "Winter"),
      "Adult",
      classification
    )
  )

# 转换为长格式并正确合并均值与标准差
time_table <- time_long %>%
  # 将户外、室内时间分别变为长格式，每条记录包含 mean 和 sd
  pivot_longer(
    cols = starts_with("Time_"),
    names_to = c("microenvironment", "stat"),
    names_pattern = "Time_(outdoor|indoor)_(log_mean|log_sd)",
    values_to = "value"
  ) %>%
  # 格式化微环境名称
  mutate(microenvironment = recode(microenvironment,
                                   outdoor = "Outdoor",
                                   indoor  = "Indoor")) %>%
  # 宽格式：log_mean 和 log_sd 各占一列
  pivot_wider(
    names_from = stat,
    values_from = value
  ) %>%
  # 添加元数据列
  mutate(
    time_dist = "lognormal",
    time_unit = "min/d"
  ) %>%
  select(classification, season, microenvironment, time_dist,
         time_log_mean = log_mean, time_log_sd = log_sd, time_unit, ref)

# ================= 保存为包内数据 =================
usethis::use_data(ir_long_table, ir_short_table, time_table, overwrite = TRUE)
