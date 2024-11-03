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
# JacobiDavidson method

## 修正方程

JacobiDavidson方法通过求解修正方程

$$\mathbf{(I - UU^H)}(\mathbf A - \tilde \lambda \mathbf I) \mathbf{(I - UU^H)t = Bt = -r}$$

逐渐改善近似的特征向量。

**注意**:

1. 实现上应保持$\mathbf r \perp \mathbf U$，否则矩阵$\mathbf B$将是刚性的并导致数值不稳定。即解向量的模过大导致正交化失效，直观上将出现伪简并现象。
2. 由于新的特征向量将在子空间$\{ \mathbf{U, t} \}$中搜索，上式中$\mathbf r$的符号不重要。

## refinedSearch()的实现

改进的Ritz搜索$^{[1]}$(Refined Ritz extraction)要求在子空间$U$中找到线性组合$\mathbf c$使得

$$\hat{\mathbf c} = \argmin_{|\mathbf c| = 1} |(AU - \lambda U)\mathbf c|$$

令$B = AU- \lambda U$, 则$$|(AU - \lambda U)\mathbf c| = \mathbf c^H B^H B \mathbf c$$

注意到$B^H B$是正定的。上式最小值即$B^H B$的最小本征值，此时$\mathbf c$为最小本征值对应的特征向量。最小化问题转化为一个小矩阵的本征值问题。

## Reference

[1] GAMM-Mitteilungen 29(2), 368-382 (2006); https://doi.org/10.1002/gamm.201490038
