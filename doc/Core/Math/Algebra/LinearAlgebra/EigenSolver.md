<!--
Copyright 2024-2025 Weibo He.

This file is part of Physica.

Permission is granted to copy, distribute and/or modify this document
under the terms of the GNU Free Documentation License, Version 1.3
or any later version published by the Free Software Foundation;
with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.

You should have received a copy of the GNU Free Documentation License
along with Physica.  If not, see <https://www.gnu.org/licenses/>.
-->
# EigenSolver - Stardard Eigenvalue Problem

## Solving eigenvectors

If we were to solve eigenvectors from real Schur decomposition, we have to solve an upper quasi-triangular equation. Solution of this equation is verbose, refer to [1] for a neat numerical implementation.

### $\lambda \in \mathbb{R}$

The equation to solve$^{[2]}$ is:

$$(T_{11} - \lambda I)w = -u$$,

which can be rewritten as

$$\left[ \begin{matrix} T & u \end{matrix} \right] \left[ \begin{matrix} w \\ 1 \end{matrix} \right] = 0$$

The last element of solution vector is known as 1. Suppose we have got the value from $w_{i + 1}$ to $w_{n}$ we want to calculate $w_i$
(1) If ${T_{ii-1}} = 0$, so we have

$$(w_i, w_{i + 1}, ..., w_n, 1) \cdot (T_{ii}, ..., T_{in}, u_i) = 0$$

$$\implies w_i = \frac{-(w_{i + 1}, ..., w_n, 1) \cdot (T_{ii + 1}, ..., T_{in}, u_i)}{T_{ii}}$$

Note that $u_{i + 1} ... u_n$ is not used to compute $w_i$, so it is efficient to store $w_i$ to $u_i$

(2) If ${T_{ii-1}} < 0$, we have

$$\left[ \begin{matrix} T_{i - 1\ i - 1} & ... & T_{i - 1\ n} & u_{i - 1} \\ T_{i\ i - 1} & ... & T_{i\ n} & u_{i} \end{matrix} \right] \left[ \begin{matrix} w_{i - 1} \\ w_i \\ ... \\ w_n \\ 1 \end{matrix} \right] = 0$$

Let

$$a = -(w_{i + 1}, ..., w_n, 1) \cdot (T_{i - 1, i + 1}, ..., T_{i - 1, n}, u_{i - 1})$$

$$b = -(w_{i + 1}, ..., w_n, 1) \cdot (T_{i, i + 1}, ..., T_{i, n}, u_{i})$$

The equation becomes to

$$\left[ \begin{matrix} T_{i - 1\ i - 1} & T_{i - 1\ i} \\ T_{i\ i - 1} & T_{i\ i} \end{matrix} \right] \left[ \begin{matrix} w_{i - 1} \\ w_i \end{matrix} \right] = \left[ \begin{matrix} a \\ b \end{matrix} \right]$$

Let's use clammer's rule to solve it. The determinate of the coefficient matrix $C$ is

$$det(C) = \left| \begin{matrix} T_{i - 1\ i - 1} & T_{i - 1\ i} \\ T_{i\ i - 1} & T_{i\ i} \end{matrix} \right| = (Re(\lambda_C)) - \lambda)^2 + Im(\lambda_C)^2$$

where $\lambda_C$ is one of eigenvalues of $C$.

### $\lambda \in \mathbb{C}$

The equation to solve is

$$ \left[ \begin{matrix} T_{11} & u_1 & u_2 & T_{14} \\ 0 & a & b & v^T_1 \\ 0 & c & d & v^T_2 \\ 0 & 0 & 0 & T_{41} \end{matrix} \right] \left[ \begin{matrix} w \\ w' \\ w'' \\ 0 \end{matrix} \right] = \lambda \left[ \begin{matrix} w \\ w' \\ w'' \\ 0 \end{matrix} \right] $$,

where value of $w'$ and $w''$ can be fetched by computing the eigenvectors of $\left[ \begin{matrix} a & b \\ c & d \end{matrix} \right]$, we have

$$w' = \frac{\frac{a - d}{2} \pm i \sqrt{|\frac{(a - d)^2}{4} + bc|}}{c}$$

$$w'' = 1$$

Now the equation becomes to

$$(T_{11} - \lambda I)w + w' u_1 + u_2 = 0$$

Complex eigenvalues corresponds to complex eigenvectors, so we have to divide the equation into two parts for the convenience of implementation. Note that $w''$ is a real number, so it is necessary to initialize vector $r$ on the right hand side of vector $c$ and the matrix becomes to a upper triangular matrix.

Let

$$w = r + ic$$

$$A = T_{11} - Re(\lambda)I$$

$$B = -Im(\lambda)I$$

The equation becomes to

$$\left \{ \begin{matrix} Ar - Bc + Re(w') u_1 + u_2 = 0 \\ Ac + Br + Im(w') u_1 = 0 \end{matrix} \right.$$

Suppose we have got values of $r_{i + 1} ... r_{n}$ and $c_{i + 1} ... c_{n}$ we want to calculate $w_i$

(1) If $A_{i i - 1} = 0$, we have

$$(A_{ii}, A_{i, i + 1}, ..., A_{in}, u_1, u_2) \cdot (r_i, r_{i + 1}, ..., r_n, Re(w'), 1) - B_{ii} c_i = 0$$

$$(A_{ii}, A_{i, i + 1}, ..., A_{in}, u_1, u_2) \cdot (c_i, c_{i + 1}, ..., c_n, Im(w'), 0) + B_{ii} r_i = 0$$

Let

$$\xi = -(A_{i, i + 1}, ..., A_{i, n}, u_1, u_2) \cdot (r_{i + 1}, ..., r_n, Re(w'), 1)$$

$$\eta = -(A_{i, i + 1}, ..., A_{i, n}, u_1, u_2) \cdot (c_{i + 1}, ..., c_n, Im(w'), 0)$$

We can get the solution:

$$ \left \{ \begin{matrix} r_i = \frac{A_{ii} \xi + B_{ii} \eta}{A_{ii}^2 + B_{ii}^2} \\ c_i = \frac{A_{ii} \eta - B_{ii} \xi}{A_{ii}^2 + B_{ii}^2} \end{matrix} \right. $$

Analogitically, it is efficient to save $r_i$ and $c_i$ to $u_1$ and $u_2$.

(2) If $A_{i i - 1} < 0$, we have

$$(A_{ii - 1}, A_{ii}, ..., A_{in}, u_1, u_2) \cdot (r_{i - 1}, r_{i}, ..., r_n, Re(w'), 1) - B_{ii} c_i = 0$$

$$(A_{ii - 1}, A_{i, i}, ..., A_{in}, u_1, u_2) \cdot (c_{i - 1}, c_{i}, ..., c_n, Im(w'), 0) + B_{ii} r_i = 0$$

$$(A_{i - 1 i - 1}, A_{i - 1 i}, ..., A_{i - 1 n}, u_1, u_2) \cdot (r_{i - 1}, r_{i}, ..., r_n, Re(w'), 1) - B_{i - 1 i - 1} c_{i - 1} = 0$$

$$(A_{i - 1 i - 1}, A_{i - 1, i}, ..., A_{i - 1 n}, u_1, u_2) \cdot (c_{i - 1}, c_{i}, ..., c_n, Im(w'), 0) + B_{i - 1 i - 1} r_{i - 1} = 0$$

Let

$$\xi_1 = -(A_{i, i + 1}, ..., A_{i, n}, u_1, u_2) \cdot (r_{i + 1}, ..., r_n, Re(w'), 1)$$

$$\xi_2 = -(A_{i, i + 1}, ..., A_{i, n}, u_1, u_2) \cdot (c_{i + 1}, ..., c_n, Im(w'), 0)$$

$$\xi_3 = -(A_{i - 1, i + 1}, ..., A_{i - 1, n}, u_1, u_2) \cdot (r_{i + 1}, ..., r_n, Re(w'), 1)$$

$$\xi_4 = -(A_{i - 1, i + 1}, ..., A_{i - 1, n}, u_1, u_2) \cdot (c_{i + 1}, ..., c_n, Im(w'), 0)$$

The equation can be written as

$$\left[ \begin{matrix} A_{i i - 1} & 0 & A_{ii} & -B_{ii} \\ 0 & A_{i i - 1} & B_{ii} & A_{ii} \\ A_{i - 1 i - 1} & -B_{i - 1 i - 1} & A_{i - 1 i} & 0 \\ B_{i - 1 i - 1} & A_{i - 1 i - 1} & 0 & A_{i - 1 i} \end{matrix} \right] \left[ \begin{matrix} r_{i - 1} \\ c_{i - 1} \\ r_i \\ c_i \end{matrix} \right] = \left[ \begin{matrix} \xi_1 \\ \xi_2 \\ \xi_3 \\ \xi_4 \end{matrix} \right]$$

## 前向自动微分

沿用[3]的符号, 对与特征值分解$\mathbf{AU = UD}$有

$$\mathbf{E} \circ \text{d}\mathbf{C} + \text{d}\mathbf{D} = \mathbf{U}^{-1} \text{d}\mathbf{AU}$$

特征值的导数:

$$\text{d}\mathbf{D} = \mathbf{I} \circ (\mathbf{U}^{-1} \text{d}\mathbf{AU})$$

特征向量的导数:

$$\text{d}\mathbf{U} = \mathbf{U}[\mathbf{F} \circ (\mathbf{U}^{-1} \text{d}\mathbf{AU})] + \mathbf{UD'}$$

其中$\mathbf{D'}$为任意对角矩阵。对于归一化的特征向量$|\mathbf{u}|^2 = 1$, 两边微分有

$$\mathbf{u} \cdot \text{d}\mathbf{u} = 0$$

即特征向量的梯度正交于特征向量自身, 可作如下正交化消除任意性

$$\text{d}\mathbf{U} \to \text{d}\mathbf{U} \cdot [\mathbf{I} - \mathbf{I} \circ \mathbf{U}^\dagger \text{d}\mathbf{U}]$$

若特征向量的模对计算没有影响, 可以简单地取$\mathbf{D'} = 0 ^{[3]}$。

## References

[1] Eigen; https://eigen.tuxfamily.org  
[2] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013  
[3] Giles, M. An extended collection of matrix derivative results for forward and reverse mode algorithmic differentiation (2008); https://people.maths.ox.ac.uk/gilesm/files/NA-08-01.pdf.
