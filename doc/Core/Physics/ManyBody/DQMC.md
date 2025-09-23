<!--
Copyright 2025 Weibo He.

This file is part of Physica.

Permission is granted to copy, distribute and/or modify this document
under the terms of the GNU Free Documentation License, Version 1.3
or any later version published by the Free Software Foundation;
with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.

You should have received a copy of the GNU Free Documentation License
along with Physica.  If not, see <https://www.gnu.org/licenses/>.
-->
# Notes on Implementation of DQMC

## 秩1更新

我们将$\Delta$矩阵从右边作用, 因此更新形式与$[1]$略有不同, 忽略自旋指标, Eq. (7.43)与Eq. (7.45)应分别修改为

$$R = 1 + (1 - G_{ii})\Delta_{ii}$$

$$G \to G - \frac{1}{R}(I - G)\Delta G$$

## Reference

[1] Gubernatis J, Kawashima N, Werner P. Quantum Monte Carlo Methods: Algorithms for Lattice Models. Cambridge University Press; 2016:194  
