<!--
Copyright 2025-2026 Weibo He.

This file is part of Physica.

Permission is granted to copy, distribute and/or modify this document
under the terms of the GNU Free Documentation License, Version 1.3
or any later version published by the Free Software Foundation;
with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.

You should have received a copy of the GNU Free Documentation License
along with Physica.  If not, see <https://www.gnu.org/licenses/>.
-->
# DiffBerry - AutoDiff approach to Berry curvature

本示例使用自动微分方法计算贝里曲率, 积分得到能带陈数。自动微分方法的优势在于只需要一个本征态, 为回避有限差分方法中相位不确定的问题提供新方法。

考虑二维晶格上的两能带陈绝缘体$^{[1]}$:

$$H(\mathbf k) = \vec{d}(\mathbf{k}) \cdot \vec{\sigma}$$

其中$\vec{d}(\mathbf{k}) = (\sin(k_x), 3\sin(k_y), 1 - \cos(k_x) - \cos(k_y))$, $\vec{\sigma} = (\sigma_x, \sigma_y, \sigma_z)$为泡利矩阵。贝里曲率

$$\Omega(\mathbf{k}) = -2\; \text{Im}[\braket{\frac{\partial \psi}{\partial k_x}|\frac{\partial \psi}{\partial k_y}}]$$

其中$\ket{\psi}$为非简并本征态。基态和激发态的陈数分别为$-1$和$1$。

## Reference

[1] Phys. Rev. B 91, 165140; https://doi.org/10.1103/PhysRevB.91.165140
