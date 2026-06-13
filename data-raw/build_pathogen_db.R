library(readr)
library(dplyr)
library(usethis)

# 手动指定列类型：所有参数列均为 double（缺失值自动转为 NA）
pathogen_db <- read_csv("data-raw/pathogen_params.csv",
                        col_types = cols(
                          pathogen_id    = col_character(),
                          pathogen_cn    = col_character(),
                          dose_response_model = col_character(),
                          alpha          = col_double(),
                          beta           = col_double(),
                          N50            = col_double(),
                          k              = col_double(),
                          dist           = col_character(),
                          alpha_log_mean = col_double(),
                          alpha_log_sd   = col_double(),
                          beta_log_mean  = col_double(),
                          beta_log_sd    = col_double(),
                          N50_log_mean   = col_double(),
                          N50_log_sd     = col_double(),
                          k_log_mean     = col_double(),
                          k_log_sd       = col_double(),
                          dose_unit      = col_character(),
                          inhalation_exposure = col_character(),
                          dermal_contact = col_character(),
                          ingestion_exposure  = col_character(),
                          response_type  = col_character(),
                          `human-to-human transmission` = col_character(),
                          ref            = col_character()
                        ),
                        na = c("", "NA", "N/A", "-")
)

# 统一模型列名和模型名称
pathogen_db <- pathogen_db %>%
  rename(model = dose_response_model) %>%
  mutate(model = recode(model,
                        "approx_beta_poisson"    = "approx_beta_poisson",
                        "exact_beta_poisson"     = "exact_beta_poisson",
                        "exponential"            = "exponential"
  ))

usethis::use_data(pathogen_db, overwrite = TRUE)
