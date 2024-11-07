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
# BKSModel

## 模板参数

**AvoidTooNear**

![](./BKSModel.png)

**图1** BKS势$^{[1]}$属于类Buckingham势，Buckingham势在$r \to 0$时将出现上图所示的不合理的吸引力，这将导致两个原子无限接近。在实践中低温下的固体SiO2不会出现原子过近的现象。但在高温的SiO2中，原子有足够的动能翻越势垒而进入左侧的吸引区，这将导致原子坐标发散。

Buckingham势为

$$V(r) = A e^{-br} - \frac{c}{r^6}$$

$$F(r) = -\frac{\partial V}{\partial r} = Ab e^{-br} - \frac{6c}{r^7}$$

实践上，低温状态下可以使用上式，为避免原子过近，高温下将Buckingham修正为

$$V_1(r) = \left\{ \begin{matrix} 2V(r_0) - V(r) \\ V(r) \end{matrix} \quad \begin{matrix} (r < r_0) \\ (r_0 \le r) \end{matrix} \right.$$

$$F_1(r) = \left\{ \begin{matrix} -F(r) \\ F(r) \end{matrix} \quad \begin{matrix} (r < r_0) \\ (r_0 \le r) \end{matrix} \right.$$

其中$r_0$应满足$F(r_0) = 0$，数值解得

$$r_0(\mathrm{O-O}) = 2.8414313142730038 \; \mathrm{Bohr}$$

$$r_0(\mathrm{Si-O}) = 2.0754014858102230 \; \mathrm{Bohr}$$

## Reference

[1] Phys. Rev. Lett. 64, 1955 (1990); https://doi.org/10.1103/PhysRevLett.64.1955
