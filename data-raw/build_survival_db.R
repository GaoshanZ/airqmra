# 从 CSV 构建活菌比例数据集
library(readr)
library(usethis)

survival_db <- read_csv("data-raw/survival_params.csv",
                        col_types = cols(
                          pathogen_id = col_character(),
                          dist        = col_character(),
                          log_mean    = col_double(),
                          log_sd      = col_double(),
                          value       = col_double(),
                          ref         = col_character()
                        ))
# 确保所有字符列标记为 UTF‑8
survival_db <- survival_db %>%
  mutate(across(where(is.character), enc2utf8))

usethis::use_data(survival_db, overwrite = TRUE)
