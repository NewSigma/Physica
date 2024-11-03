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
# FokkerPlanck

本示例演示如何使用有限元方法求解有阻尼的谐振子的Fokker-Planck方程$^{[1]}$:

$$\frac{\partial W(x, p, t)}{\partial t} = [-\frac{p}{m} \frac{\partial}{\partial x} + \frac{\partial}{\partial p} (x + \gamma p) + \frac{\partial^2}{\partial p^2} D] W(x, p, t)$$

其中谐振子质量为$m$。$x$和$p$分别为谐振子的位置和动量。$\gamma$为阻尼系数，$D$为扩散系数。$W$为系统的概率密度分布。

![](FokkerPlanck.gif)

**图1** 概率密度$W(x, p, t)$随时间演化。x轴表示位置，y轴表示动量。可视化使用Paraview$^{[2]}$完成

## Reference

[1] Nonlinear Dyn. 4, 357–372 (1993); https://doi.org/10.1007/BF00120671
[2] Paraview; https://www.paraview.org/
