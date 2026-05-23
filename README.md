# airqmra: 吸入暴露定量微生物风险评估 R 包

`airqmra` 是一个用于**大气环境生物气溶胶感染风险评估**的专业 R 包。
它集成了中国人群暴露参数数据库、多病原剂量‑反应模型、一维与二维蒙特卡洛模拟、年风险与疾病负担计算，以及多层次敏感性分析工具（Spearman、PRCC、Morris、Sobol'）。
包内同时提供了交互式 Shiny 网页应用和符合团体标准（T/CSES XXXX-XXXX）的自动化报告生成。

## 安装

```r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("GaoshanZ/airqmra")
```

## 快速示例

```r
library(airqmra)

# 创建暴露场景：成人室外吸入嗜肺军团菌
scen <- new_scenario(
  pathogen_id         = "legionella pneumophila",
  exposure_group      = "Adult",
  microenvironment    = "Outdoor",
  population_deposition = "Male",
  use_uncertainty     = FALSE
)

# 一维蒙特卡洛模拟（仅变异性）
single_risk <- run_qmra_1d(scen, n_iter = 5000)
annual      <- annual_risk(single_risk, frequency = 20)

# 二维蒙特卡洛（变异性 + 参数不确定性）
scen_unc <- new_scenario(
  pathogen_id         = "legionella pneumophila",
  exposure_group      = "Adult",
  microenvironment    = "Outdoor",
  population_deposition = "Male",
  use_uncertainty     = TRUE
)
res2d <- run_qmra_2d(scen_unc, n_var = 2000, n_unc = 500, frequency = 20)

# 敏感性分析
sens_spearman <- sensitivity_1d(scen, n_iter = 3000)
tornado_plot(sens_spearman, title = "Spearman sensitivity")

prcc_res <- sensitivity_prcc(scen, n_iter = 3000, nboot = 200)
tornado_plot(prcc_res, title = "PRCC sensitivity")

mor <- sensitivity_morris(scen, n_levels = 8, n_trajectories = 30)
tornado_plot(mor, title = "Morris sensitivity (mu*)", log_scale = TRUE)

sobol_res <- sobol_sensitivity(scen_unc, n = 1000)
plot(sobol_res)

# 生成符合团标的风险评估报告
generate_risk_report(scen, list(single_risk = single_risk,
                                annual_risk  = annual,
                                sensitivity_1d = sens_spearman),
                     language = "zh")
```

## 主要功能

### 📦 内置数据库
- 15+ 种常见病原体的剂量‑反应参数（指数、精确 β‑泊松、近似 β‑泊松），包含不确定性分布（对数正态）。
- 中国人群暴露参数：长期/短期呼吸速率、室内外时间分配、肺部沉积率、气溶胶存活率、活动情景模板。
- 所有参数支持外部 CSV 驱动，用户可自由扩展。

### 🔬 剂量‑反应模型
- 指数模型、精确 β‑泊松模型（数值安全保护）、基于 N₅₀ 的近似 β‑泊松模型。
- **模型拟合模块**：支持从原始实验数据（剂量‑反应表）进行最大似然拟合、bootstrap 参数分位数、模型比较（似然比检验）。

### 🎲 蒙特卡洛模拟
- **一维模拟**：传播暴露变异性（浓度、吸入速率、时间、沉积率、存活率）。
- **二维模拟**：分离参数不确定性（k, α, β, N₅₀）与人群变异性，输出风险均值的置信区间及综合风险分布。
- 年感染风险计算（支持固定频率或频率分布）。

### 📊 敏感性分析
- **Spearman 秩相关**（1D 变异性 / 2D 参数不确定性 / 全局混合）。
- **PRCC（偏秩相关系数）**：控制其他变量后的单调影响，含 bootstrap 置信区间。
- **Morris 基本效应法**：自包含实现，筛选重要变量及交互效应。
- **Sobol' 方差分解**（基于 sobol2007 算法）：一阶与总效应指数，完整主效应与交互效应量化。
- 飓风图 (`tornado_plot`) 及标准散点图可视化。

### 🌐 数据交互与验证
- 浓度数据导入：支持 CSV/TXT，自动拟合分布（对数正态、Gamma、Weibull），并提供数据诊断。
- 场景检查：`check_scenario()` 自动评估参数合理性。
- 灵活的 `update_scenario()` 方法，便于批量情景分析。

### 📄 报告生成
- 一键生成中英文 HTML/Word/PDF 风险评估报告，严格遵循《大气环境生物气溶胶健康风险评估技术方法指南》(T/CSES XXXX-XXXX) 框架。
- 报告内容自适应：根据用户提供的分析结果，自动展示 Spearman、PRCC、Morris、Sobol' 及二维不确定性分析。

### 🖥️ Shiny 交互界面
- 内置完整的网页应用 (`run_shiny()`)，支持场景设置、一维/二维模拟、多敏感性分析、情景比较和报告下载。

## 支持的部分病原体

| 病原体 | 模型 | 点估计参数 |
|--------|------|-------------|
| 嗜肺军团菌 | 指数 | k = 4.99×10⁻⁵ |
| 金黄色葡萄球菌 | 指数 | k = 7.64×10⁻⁸ |
| 铜绿假单胞菌 | 近似 β‑泊松 | α = 0.19, N₅₀ = 1.85×10⁴ |
| 烟曲霉菌 | 精确 β‑泊松 | α = 0.09, β = 14 |
| 流感病毒 | 近似 β‑泊松 | α = 0.581, N₅₀ = 9.45×10⁵ |

（完整列表请参见包内 `pathogen_db` 数据集）

## 依赖
- `dplyr`, `tidyr`, `readr` – 数据处理
- `hypergeo` – 精确 β‑泊松模型
- `sensitivity` – Sobol' 分析
- `shiny`, `shinythemes` – 交互界面
- `rmarkdown`, `knitr` – 报告生成
- `ggplot2`, `scales` – 可视化
- R ≥ 4.0

## 引用

若使用本包进行研究，请引用：

```bibtex
@software{airqmra,
  author = {Gaoshan Zhang},
  title  = {airqmra: Inhalation Quantitative Microbial Risk Assessment in R},
  year   = {2025},
  url    = {https://github.com/GaoshanZ/airqmra},
  note   = {R package version 0.0.0.9000}
}
```

暴露参数及病原体数据主要来源：
- Duan, X.L. et al. (2014, 2016) *Highlights of the Chinese Exposure Factors Handbook*.
- [QMRA Wiki](https://qmrawiki.msu.edu) 及公开文献。

## 贡献与反馈

欢迎通过 [GitHub Issues](https://github.com/GaoshanZ/airqmra/issues) 提交问题或建议。
作者邮箱：2024029025@chd.edu.cn

---

**让风险分析更透明，从黑箱走向科学可重复！**
