<!--
Copyright 2024 Weibo He.

This file is part of Physica.

Permission is granted to copy, distribute and/or modify this document
under the terms of the GNU Free Documentation License, Version 1.3
or any later version published by the Free Software Foundation;
with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.

You should have received a copy of the GNU Free Documentation License
along with Physica.  If not, see <https://www.gnu.org/licenses/>.
-->
# ProbDistribution

ProbDistribution使用采样语义计算一维随机变量的概率分布函数, 考虑随机变量

$$y = \ln\left( 1 + \sqrt{1-e^{-1/n}} \right) \sum_{i = 1}^n s_i$$

其中独立随机变量$s_i = \pm 1$满足$p = 1/2$的二项分布

![](./ProbDistribution.png)

**图1** 蓝色: 随机变量$y$的概率密度分布; 绿色: 标准正态分布

随机变量$y$的各阶矩与标准正态分布的随机变量相等, 且矩母函数存在, 故$n \to \infty$时$y \sim N(0, 1)$
