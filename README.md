# airqmra: 吸入暴露定量微生物风险评估 R 包

`airqmra` 是一个专门用于**吸入暴露病原微生物定量风险评估 (QMRA)** 的 R 包。  
它集成了中国人群暴露参数数据库、多种病原体的剂量‑反应模型、一维与二维蒙特卡洛模拟、年风险计算以及全面的敏感性分析工具。  
本包力求在可重复性、透明度和灵活性上超越商业软件 **Crystal Ball**，为科研和公共卫生决策提供可靠的开源工具。

## 安装

```r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("GaoshanZ/airqmra")
快速示例
r
library(airqmra)

# 创建一个暴露场景：成人室外吸入嗜肺军团菌
scen <- new_scenario(
  pathogen_id         = "legionella pneumophila",
  exposure_group      = "Adult",
  microenvironment    = "Outdoor",
  population_deposition = "Male",
  use_uncertainty     = FALSE      # 点估计模式
)

# 一维蒙特卡洛模拟（仅变异性）
single_risk <- run_qmra_1d(scen, n_iter = 5000)
summary(single_risk)

# 计算年风险（每年暴露 20 次）
annual <- annual_risk(single_risk, frequency = 20)
summary(annual)

# 二维蒙特卡洛（变异性 + 参数不确定性）
scen_unc <- new_scenario(
  pathogen_id         = "legionella pneumophila",
  exposure_group      = "Adult",
  microenvironment    = "Outdoor",
  population_deposition = "Male",
  use_uncertainty     = TRUE
)
res2d <- run_qmra_2d(scen_unc, n_var = 2000, n_unc = 500, frequency = 20)
summary(res2d)
plot(res2d)

# 敏感性分析（Spearman 飓风图）
sens1d <- sensitivity_1d(scen, n_iter = 3000)
tornado_plot(sens1d, title = "1D Sensitivity")

# Sobol' 方差分解（全局敏感性）
sobol_res <- sobol_sensitivity(scen_unc, n = 500)
plot(sobol_res)
主要功能
内置数据库

14 种常见病原体的剂量‑反应参数（指数、精确 β‑泊松、近似 β‑泊松三种模型），并包含参数不确定性信息（对数正态分布）。

中国人群暴露参数：长期/短期呼吸速率、室内外时间分配、肺部沉积率、气溶胶存活率。

活动情景模板（上班族、户外工人、学生、老人等）。

剂量‑反应模型

指数模型

精确 β‑泊松模型（数值安全保护）

基于 N₅₀ 的近似 β‑泊松模型

蒙特卡洛模拟

一维模拟：传播暴露变异性（浓度、吸入速率、时间、沉积率、存活率）。

二维模拟：分离参数不确定性（k, α, β, N₅₀）与人群变异性，输出风险均值的置信区间及综合风险分布。

年感染风险计算（支持固定频率或频率分布）。

敏感性分析

Spearman 秩相关（1D 变异性 / 2D 参数不确定性 / 全局混合）。

Sobol' 方差分解（完整主效应与交互效应，一阶及总效应指数）。

飓风图可视化 (tornado_plot)。

灵活的场景搭建

通过关键词（人群、季节、微环境、病原体 ID）自动加载默认参数。

所有输入均可由用户自定义分布函数。

支持的部分病原体
病原体	模型	点估计参数
嗜肺军团菌	指数	k = 4.99×10⁻⁵
金黄色葡萄球菌	指数	k = 7.64×10⁻⁸
铜绿假单胞菌	近似 β‑泊松	α = 0.19, N₅₀ = 1.85×10⁴
烟曲霉菌	精确 β‑泊松	α = 0.09, β = 14
流感病毒	近似 β‑泊松	α = 0.581, N₅₀ = 9.45×10⁵
（完整列表请参见包内 pathogen_db 数据集）

依赖
dplyr, tidyr, readr – 数据处理

hypergeo – 精确 β‑泊松模型

sensitivity – Sobol' 分析

ggplot2（可选，用于部分绘图）

R ≥ 4.0

引用
若使用本包进行研究，请引用：

bibtex
@software{airqmra,
  author = {Gaoshan Zhang},
  title  = {airqmra: Inhalation Quantitative Microbial Risk Assessment in R},
  year   = {2025},
  url    = {https://github.com/GaoshanZ/airqmra},
  note   = {R package version 0.0.0.9000}
}
暴露参数及病原体数据主要来源：

Duan, X.L. et al. (2014, 2016) Highlights of the Chinese Exposure Factors Handbook.

QMRA Wiki 及公开文献。

贡献与反馈
欢迎通过 GitHub Issues 提交问题或建议。
作者邮箱：2024029025@chd.edu.cn

让风险分析更透明，从黑箱走向科学可重复！