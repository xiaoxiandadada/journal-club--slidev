# journal club --slidev

## slides 制作工具
- 基于 Slidev 制作了这个演示文档，Slidev是一个专为开发者打造的演示文稿工具，基于 Vite 构建，并使用 Vue 3 作为核心框架，可以使用熟悉的 Markdown 语法来制作炫酷的 PPT，它还支持 HTML 和 Vue 组件，能够呈现像素级完美的布局。它的GitHub 仓库 [https://github.com/slidevjs/slidev](https://github.com/slidevjs/slidev)

## 文章主要内容
- 讲解的这篇文章《A scalable variational inference approach for increased mixed-model association power》介绍了一种名为Quickdraws的可扩展变分推断方法，旨在提高混合模型关联分析的统计能力，同时保持计算效率。
### 方法介绍
Quickdraws方法通过以下方式提高关联分析能力：
- **贝叶斯回归模型**：采用尖峰和板先验（spike-and-slab prior）对变异效应进行建模，这种先验能够更好地捕捉非多基因性状的结构。
- **随机变分推断**：使用随机变分推断（stochastic variational inference）来高效训练模型，结合一阶优化器和图形处理单元（GPU）矩阵运算，实现线性扩展的计算复杂度。
- **模型拟合与测试**：在模型拟合阶段，Quickdraws首先估计遗传和环境方差组分，然后使用留一染色体外（LOCO）方案预测表型。在测试阶段，使用估计的遗传效应计算线性或逻辑混合模型的得分检验统计量。

### 实验结果
- **模拟数据**：在模拟数据中，Quickdraws在定量性状上的关联能力与BOLT-LMM相当，但计算需求显著减少。对于二元性状，其关联能力高于SAIGE、REGENIE和FastGWA-GLMM。
- **真实数据**：在分析英国生物库（UK Biobank）中约40.5万名个体的79个定量性状和50个疾病性状时，Quickdraws比REGENIE和FastGWA分别多检测到4.97%和22.71%的关联，对于疾病性状则多检测到3.25%和7.07%的关联。

### 贡献和影响
- **提高GWAS效率**：Quickdraws在保持计算效率的同时提高了关联分析能力，为大规模GWAS提供了更有效的工具。
- **提供新分析工具**：为遗传学家提供了一种新的分析工具，能够更好地处理复杂的遗传数据，发现更多的遗传关联。
- **推动机器学习在遗传学中的应用**：展示了机器学习技术在遗传学研究中的潜力，为未来的研究提供了新的思路和方法。
