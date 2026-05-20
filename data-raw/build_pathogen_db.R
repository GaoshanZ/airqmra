# 读取并处理病原体参数 CSV，生成 pathogen_db.rda
library(readr)
library(dplyr)
library(usethis)

pathogen_db <- read_csv("data-raw/pathogen_params.csv",
                        col_types = cols(.default = col_character()))

# 数值列转换
numeric_cols <- c("alpha", "beta", "N50", "k")
pathogen_db <- pathogen_db %>%
  mutate(across(any_of(numeric_cols), as.numeric))

# 统一模型名称（确保与 dose_response 匹配）
pathogen_db <- pathogen_db %>%
  mutate(model = recode(model, "beta_poisson" = "exact_beta_poisson"))

# 由 N50 推算缺失的 beta（仅用于 β‑泊松模型且 alpha、N50 已知）
pathogen_db <- pathogen_db %>%
  mutate(
    beta = ifelse(
      (model == "exact_beta_poisson") & is.na(beta) & !is.na(alpha) & !is.na(N50),
      N50 / (2^(1/alpha) - 1),
      beta
    )
  )

# 保存为包数据
usethis::use_data(pathogen_db, overwrite = TRUE)
