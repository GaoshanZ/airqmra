# 病原体剂量‑反应参数数据库（英文主键版）
library(tibble)
library(dplyr)

# 原始数据逐条录入
pathogen_db_raw <- tibble::tribble(
  ~pathogen_cn,         ~pathogen_id,                         ~model,         ~alpha, ~beta, ~N50,     ~k,       ~primer_F,                          ~primer_R,
  "大肠杆菌O157:H7",   "escherichia coli o157:h7",           "beta_poisson",  0.248, 48.8,  NA,       NA,       NA,                                 NA,
  "烟曲霉菌",          "aspergillus fumigatus",              "beta_poisson",  1.1,   20,    NA,       NA,       NA,                                 NA,
  "产肠毒素大肠杆菌",  "enterotoxigenic escherichia coli",   "beta_poisson",  0.178, NA,    8.6e7,    NA,       "TCCTGAAAGCATGAATAGTAGCAATTAC",    "TTAATAGCACCCGGTACAAGCA",
  "福氏志贺氏菌",      "shigella flexneri",                  "beta_poisson",  0.265, NA,    1.48e3,   NA,       "ACCATGCTCGCAGAGAAACT",            "TACGCTTCAGTACAGCATGC",
  "空肠弯曲杆菌",      "campylobacter jejuni",               "beta_poisson",  0.151, NA,    1.69e3,   NA,       "CTGGTGGTTTTGAAGCAAAGATT",         "CAATACCAGTGTCTAAAGTGCGTTTAT",
  "金黄色葡萄球菌",    "staphylococcus aureus",              "exponential",   NA,    NA,    NA,       7.64e-8,  NA,                                 NA,
  "金黄色葡萄球菌",    "staphylococcus aureus",              "exponential",   NA,    NA,    NA,       8.05e-8,  NA,                                 NA,
  "志贺氏杆菌",        "shigella spp.",                      "beta_poisson",  0.265, NA,    1.48e3,   NA,       NA,                                 NA,
  "绿脓杆菌",          "pseudomonas aeruginosa",             "beta_poisson",  0.19,  NA,    1.85e4,   NA,       "AGATGGGCGACGAGTTGCTG",            "GACTGATGACCGTGGGCTTGAT",
  "弯曲杆菌属",        "campylobacter spp.",                 "beta_poisson",  0.145, NA,    8.96e2,   NA,       "CACGTGCTACAATGGCATAT",            "GGCTTCATGCTCTCGAGTT",
  "大肠杆菌O157:H7",   "escherichia coli o157:h7",           "beta_poisson",  0.209, NA,    1.12e3,   NA,       "TTTCACACTTATTGGATGGTCTCAA",       "CGATGAGTTTATCTGCAAGGTGAT",
  "金黄色葡萄球菌",    "staphylococcus aureus",              "exponential",   NA,    NA,    NA,       8.05e-8,  "TCGCTTGCTATGATTGTGG",             "ACATACGCCAATGTTCTACC",
  "肠球菌",            "enterococcus spp.",                   "exponential",   NA,    NA,    NA,       1.31e-7,  "AGAAATTCCAAACGAACTTG",            "CAGTGCTCTACCTCCATCATT",
  "铜绿假单胞菌",      "pseudomonas aeruginosa",             "beta_poisson",  0.19,  NA,    1.85e4,   NA,       "AGATGGGCGACGAGTTGCTG",            "GACTGATGACCGTGGGCTTGAT",
  "非伤寒沙门氏菌",    "nontyphoidal salmonella",            "beta_poisson",  0.3126,NA,    23618.77, NA,       "CTCACCAGGAGATTACAACATGG",         "AGCTCAGACCAAAAGTGACCATC",
  "大肠杆菌",          "escherichia coli",                    "beta_poisson",  0.248, NA,    749.51,   NA,       "GGTAGAGCACTGTTTTGGCA",            "TGTCTCCCGTGATAACTTTCTC",
  "空肠弯曲杆菌",      "campylobacter jejuni",               "beta_poisson",  0.145, 7.59,  NA,       NA,       NA,                                 NA,
  "嗜肺军团菌",        "legionella pneumophila",              "exponential",   NA,    NA,    NA,       0.06,     NA,                                 NA,
  "非伤寒沙门氏菌",    "nontyphoidal salmonella",            "beta_poisson",  0.3126,2884,  NA,       NA,       "CTCACCAGGAGATTACAACATGG",         "AGCTCAGACCAAAAGTGACCATC",
  "金黄色葡萄球菌",    "staphylococcus aureus",              "exponential",   NA,    NA,    NA,       6.46e-8,  NA,                                 NA,
  "金黄色葡萄球菌",    "staphylococcus aureus",              "exponential",   NA,    NA,    NA,       1.00e-7,  NA,                                 NA,
  "大肠杆菌",          "escherichia coli",                    "exponential",   NA,    NA,    NA,       1.22e-8,  NA,                                 NA,
  "大肠杆菌",          "escherichia coli",                    "beta_poisson",  0.155, 2.11e6,NA,       NA,       NA,                                 NA
)

# 处理数据：计算缺失的 beta（由 N50），统一模型名称
pathogen_db <- pathogen_db_raw %>%
  mutate(
    beta = ifelse(
      (model == "beta_poisson" | model == "exact_beta_poisson") & is.na(beta),
      N50 / (2^(1/alpha) - 1),
      beta
    ),
    model = recode(model,
                   "beta_poisson" = "exact_beta_poisson")
  ) %>%
  select(pathogen_id, pathogen_cn, model, alpha, beta, N50, k, primer_F, primer_R) %>%
  arrange(pathogen_id, model)

# 保存为包内数据
usethis::use_data(pathogen_db, overwrite = TRUE)
