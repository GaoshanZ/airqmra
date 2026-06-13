# 从 CSV 构建活动情景模板数据集
library(readr)
library(usethis)

activity_templates <- read_csv(
  "data-raw/activity_templates.csv",
  col_types = cols(
    template_id          = col_character(),
    Sleeping             = col_double(),
    Sedentary            = col_double(),
    Light_intensity      = col_double(),
    Moderate_intensity   = col_double(),
    High_intensity       = col_double(),
    Very_High_intensity  = col_double(),
    time_unit            = col_character()
  )
)
activity_templates <- activity_templates %>%
  mutate(across(where(is.character), enc2utf8))

# 保存为包内数据
usethis::use_data(activity_templates, overwrite = TRUE)
