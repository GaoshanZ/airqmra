# 从 CSV 构建沉积分数数据集
library(readr)
library(dplyr)
library(usethis)

deposition_db <- read_csv("data-raw/deposition_params.csv",
                          col_types = cols(
                            Population      = col_character(),
                            site            = col_character(),
                            dist            = col_character(),
                            log_mean        = col_double(),
                            log_sd          = col_double(),
                            `Deposition Site` = col_character()
                          )) %>%
  rename(
    population      = Population,
    deposition_site = `Deposition Site`
  ) %>%
  select(population, site, dist, log_mean, log_sd, deposition_site)

# 确保所有字符列标记为 UTF‑8
deposition_db <- deposition_db %>%
  mutate(across(where(is.character), enc2utf8))

usethis::use_data(deposition_db, overwrite = TRUE)
