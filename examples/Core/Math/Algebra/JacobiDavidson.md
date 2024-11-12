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
# JacobiDavidson

## Introduction

众所周知，大规模稀疏本征值问题(Sparse Eigenvalue Problem, SEP)在科学和工程领域有极为重要和广泛的应用。求解SEP的数值算法主要包括幂法$^{[1]}$, Lanzcos$^{[1]}$, Arnoldi$^{[2]}$, Davidson$^{[3]}$, Jacobi-Davidson$^{[4]}$, FEAST$^{[5]}$, LOBPCG$^{[6]}$。其中幂法只能求解极端特征值，适用范围有限。Lanczos算法可被视作Arnoldi算法在厄米情形下的特例。Arnoldi算法是目前主流科学计算软件Matlab$^{[7]}$, Mathematica$^{[8]}$所使用的算法。基于幂法改进而来的Lanczos算法和Arnoldi算法在处理内部本征值问题时将遇到困难。Davidson算法是求解对角占优矩阵的快速算法，但对于非对角占优的矩阵将会遇到难以收敛的问题。Jacobi-Davidson算法将Jacobi正交补修正技术引入Davidson算法，在保持Davidson算法高效的同时解决了Davidson算法只适用于对角占优矩阵的问题$^{[4]}$。FEAST算法用于求解任意给定区间中所有的本征值，但速度慢于Jacobi-Davidson算法。Mathematica中实现了FEAST算法。LOBPCG算法不适用于非对称矩阵。LOBPCG和JDCG(Jacobi-Davidson针对对称矩阵的变式)的作者均称自己的算法优于对方$^{[6, 10]}$, 这不免令人感到困惑。

Jacobi-Davidson算法具有以下优点:
- 速度快：Jacobi-Davidson算法收敛速度快于Davidson, Lanczos和FEAST$^{[4, 9]}$
- 适用范围广：Jacobi-Davidson算法可适用于极端特征值和内部特征值。不同于Davidson和LOBPCG，Jacobi-Davidson对矩阵类型没有要求
- 灵活性强：Jacobi-Davidson算法可分解为若干独立的模块，各种模块可自由组合以适应问题

## Main

考虑本征值问题

$$\mathbf{Ax} = \lambda \mathbf{x}$$

其中特征向量$\mathbf{x}$和特征值$\lambda$待求解。我们可以随机生成一个向量$\mathbf{x_0}$并带入上式，此时等式必不成立，有残差

$$\mathbf{r}_0(\lambda) = \mathbf{Ax}_0 - \lambda \mathbf{x}_0 \neq \mathbf{0}$$

此时$\lambda$仍是未知的。我们可以做得比随机猜测更好，即使用矩阵的Rayleigh商，它是最小二乘意义下对特征值的最优估计

$$\theta_0 = \argmin_\lambda |\mathbf{r}_0(\lambda)| = \frac{\mathbf{x}_0^T \mathbf{Ax}}{\mathbf{x}_0^T \mathbf{x}}$$
 
显然$\mathbf{r}_0(\theta_0)$与$\mathbf{x_0}$正交。随机猜测的$\mathbf{x_0}$并不能让我们满意，下面我们采用Jacobi正交补修正逐步改善解的精度。

我们记精确解可以通过修正得到$\mathbf{x} = \mathbf{x}_k + \delta \mathbf{x}_k$, $\lambda = \theta_k + \delta \theta_k$, 其中$\mathbf{x}_k, \theta_k$是第$k$次迭代得到的特征向量和特征值, 将上式带入本征值问题方程

$$\mathbf{A} (\mathbf{x}_k + \delta \mathbf{x}_k) = (\theta_k + \delta \theta_k) (\mathbf{x}_k + \delta \mathbf{x}_k)$$
 
整理并忽略二阶小量$\delta \theta_k \cdot \delta \mathbf{x}_k$有

$$(\mathbf{A} - \theta I) \delta \mathbf{x}_k - \mathbf{x}_k \delta \theta_k = -\mathbf{r}_k(\theta_k)$$

为了求解该方程，我们可以在两边乘以一个投影矩阵$I - \mathbf{x}_k^T \mathbf{x}_k$，第二项投影后结果为$0$

$$(I - \mathbf{x}_k^T \mathbf{x}_k) (\mathbf{A} - \theta I) \delta \mathbf{x}_k = -\mathbf{r}_k(\theta_k)$$
 

显然沿着$\mathbf{x}_k$方向对修正近似解没有帮助，我们感兴趣的$\delta \mathbf{x}_k$必定与$\mathbf{x}_k$正交，将解空间限制在$\mathbf{x}_k$的正交补空间:

$$(I - \mathbf{x}_k^T \mathbf{x}_k) (\mathbf{A} - \theta I) (I - \mathbf{x}_k^T \mathbf{x}_k) \delta \mathbf{x}_k = -\mathbf{r}_k(\theta_k)$$

此即Jacobi-Davidson修正方程，修正可以如下进行直至结果收敛

$$\mathbf{x}_{k + 1} = \mathbf{x}_k + \delta \mathbf{x}_k$$

$$\theta_{k + 1} = \argmin_\lambda |\mathbf{r}_{k + 1}(\lambda)|$$
 
物理上，我们通常感兴趣的是基态和最低的几个激发态，即厄米矩阵$\mathbf{A}$最小的若干本征值。Y. Notay针对这种情况提出了Jacobi-Davidson算法的一种变式，即使用共轭梯度算法(CG)求解Jacobi-Davidson修正方程。如果要求解的不是厄米矩阵或感兴趣的不是极端特征值，只需将线性求解器模块切换为MINRES、GMRES等更通用的线性求解器。这种模块的可替换性是Jacobi-Davidson算法的优点之一。

## Example

使用JDCG求解5000维厄米矩阵的前4个最小的本征值, 结果为

(-1062.865585527699, -456.5290883347747, -291.9029897722832, -216.657138458762)
