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
# Schur

## 关于splitOffTwoRows()的实现

Given a 2x2 matrix:
$$ A = \left[ \begin{matrix} a & b \\ c & d \end{matrix} \right] ,$$
whose eigenvalues are $ \frac{a + d}{2}  \pm \sqrt{\frac{(a - d)^2}{4} + bc}$

Now we want a orthogonal matrix $U$, such that

$$UAU^T = M$$

where $M$ is a upper triangular matrix.

Suppose $U$ is a givens matrix, that is:

$$ U = \left[ \begin{matrix} \cos{\theta} & \sin{\theta} \\ -\sin{\theta} & \cos{\theta} \end{matrix} \right] $$

So the (2, 1) element of $M$ is

$$ M(2, 1) = c \cos^2{\theta} - b \sin^2{\theta} + (d - a) \sin{\theta} \cos{\theta} $$

Let $M(2, 1) = 0$, note that $\sin{\theta}$ cannot be 0, or matrix $U$ will be identity, we have

$$ c \cot^2{\theta} + (d - a) \cot{\theta} - b = 0 $$

Solve the equation, we have
$$\cot{\theta} = \frac{\frac{a - d}{2}  \pm \sqrt{\frac{(a - d)^2}{4} + bc}}{c}$$

which is similiar to the expression of eigenvalues.

## 关于francisQR()的实现

考虑3x3矩阵

$$ A = \left[ \begin{matrix} a & \epsilon & 0 \\ \epsilon & a & \epsilon \\ 0 & \epsilon & a \end{matrix} \right]$$

其中$a \gg \epsilon$，需要使用QR算法将$A_{21}$化为0，令

$$s = A_{22} + A_{33} = 2a \qquad t_1 = A_{22} * A_{33} = a^2 \qquad t_2 = A_{23} * A_{32} = \epsilon^2$$

使用下式计算约化向量$c$的第一个元素

$$c_0 = (A_{00} - s) * A_{00} + t_1 + (A_{01} * A_{10} - t_2)$$

$$= (a - 2a) * a + a^2 + (\epsilon^2 - \epsilon^2)$$

注意这种括号的配置避免了$a$和$\epsilon$的直接加减运算，可以减少浮点误差的影响从而改进算法的稳定性。
