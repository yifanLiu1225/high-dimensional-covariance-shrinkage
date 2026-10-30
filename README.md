# High-Dimensional Covariance Shrinkage

## Project core

**高维协方差收缩：样本协方差的风险、谱失真与稳定化。**

本项目以理论探索为主、模拟与小规模应用为辅，研究当维数 p 与样本量 n 同阶时，样本协方差矩阵为什么会出现谱失真、不可逆和数值不稳定，以及线性收缩估计如何通过偏差—方差权衡改善估计风险。

## Central question

当 p/n 不再接近零时，普通样本协方差矩阵 S 仍然可靠吗？对于

\[
\widehat\Sigma_\lambda=(1-\lambda)S+\lambda T,
\]

如何推导最优收缩强度，解释它对风险、特征值、条件数和下游二次型的影响？

## Main theory

- Wishart 分布与样本协方差的无偏性；
- 协方差矩阵估计的 Frobenius 风险；
- 线性收缩估计的偏差—方差分解；
- 最优收缩强度的推导；
- 目标矩阵选择与正定性；
- 高维样本协方差的特征值扩散；
- Marchenko–Pastur 型谱边界的直观解释；
- 收缩估计对马氏距离、Hotelling T² 和线性判别的影响。

## Project structure

R/                 可复用 R 函数与模拟脚本
report/            理论报告与实验报告
data/raw/          原始数据或固定快照
data/processed/    处理后的分析数据
figures/           图表
tables/            结果表
references.bib     BibTeX 文献
project_plan.Rmd   详细项目大纲

## Scope

本项目不以 LLM benchmark 为主线，也不研究去偏 LASSO、Graphical LASSO、秩回归或符号秩方法。重点是协方差矩阵本身的估计风险和谱稳定性。精度矩阵、稀疏图结构和稀疏回归只在讨论中作为边界与对照，不作为主要方法。

