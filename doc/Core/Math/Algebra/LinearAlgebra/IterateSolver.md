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
# 线性系统迭代求解器

对于线性系统

$$\mathbf{Ax = b}$$

迭代近似解$\tilde{\mathbf x}$的误差为$\Delta \mathbf x = \mathbf x - \tilde{\mathbf x}$, 要求误差远小于近似解

$$\Delta \mathbf x_a \ll \tilde{\mathbf x}_a$$

其中$a$为元素指标。上式必要条件为

$$|\mathbf{A} \Delta \mathbf x| \ll |\mathbf{A} \tilde{\mathbf x}| \approx |\mathbf b|$$

利用残差$\mathbf r = \mathbf{b - A} \tilde{\mathbf x} = \mathbf{A} \Delta \mathbf x$，上式可写作

$$|\mathbf r| < \epsilon |\mathbf b|$$

其中$\epsilon$为一小量，可用来控制解的相对误差。
