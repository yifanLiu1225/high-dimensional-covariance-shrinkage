# Beyond the Leaderboard

## Project core

**超越排行榜：用多元统计识别大语言模型的潜在能力结构、评测冗余与排名稳定性。**

本项目把公开的 LLM benchmark 结果看作一个多元测量问题：模型是观测，benchmark 是变量，PCA 和因子模型用于提取潜在能力，Mahalanobis 距离用于识别能力画像异常，模拟和 bootstrap 用于验证理论与排名稳定性。

## Main data

- Hugging Face OpenEvals/leaderboard-data
- benchmark 子集和评分扰动敏感性

数据快照、版本号和处理规则必须记录，不能直接依赖实时 leaderboard。

## Planned modules

- 潜变量测量模型和理论推导
- 相关 benchmark 的信息冗余
- PCA、因子分析和模型能力画像
- 经典与稳健 Mahalanobis 距离
- PCA \(T^2\) 与 \(Q\) 诊断
- 模拟污染实验
- bootstrap 排名稳定性
- benchmark 子集和评分扰动敏感性

详细项目大纲、研究假设、数据规则、图表清单和四人分工见：

- project_plan.Rmd

## Proposed folder structure

R/                 R scripts
report/            R Markdown reports
data/raw/          frozen source files
data/processed/    cleaned analysis data
figures/           figures for the report and presentation
tables/            reproducible result tables
references/        citation files and notes

项目暂不上传 GitHub。正式分析开始前，先固定数据快照并确认每个组员的模块边界。
