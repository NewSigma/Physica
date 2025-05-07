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
# RPMD

## Overview

一般的分子动力学模拟，由四个核心概念构成:

    运动学模型(KineticModel): 运动学方程、积分算法、边界条件、构型约束条件
    相互作用模型(ForceModel): 力场
    热浴(Thermostat): 控温算法, NVT、NPT
    压浴(Barostat): 控压算法, NPT only

四者由class RPMD有机集成。以一维体系为例, RPMD哈密顿量为

$$H_n = \sum_i^N \sum_j^n [\frac{[p_i^{(j)}]^2}{2m_i} + \frac{1}{2} m_i \omega_n^2 (q_i^{(j)} - q_i^{(j + 1)})^2] + \sum_j^n V(\mathbf{q}^{(j)})$$

其中$n$为NumReplica, 当$n = 1$时, 哈密顿量退化为经典分子动力学(MD)的哈密顿量。class RPMD致力于将MD和RPMD在一个统一的框架下实现。

正则系综固定边界条件下的配分函数为

$$Q = \frac{1}{(2\pi \hbar)^{Nn}} \int e^{-\beta_n H_n(\mathbf{p}, \mathbf{q})} \mathrm{d} \mathbf{p} \mathrm{d} \mathbf{q}$$

正则系综周期性边界条件下的配分函数为

$$Q = \frac{1}{(2\pi \hbar)^{Nn}} \int e^{-\beta_n H_n(\mathbf{p}, \mathbf{q})} \delta(\sum_{ij} p_i^{(j)}) \mathrm{d} \mathbf{p} \mathrm{d} \mathbf{q}$$

## 动能计算

(1) RPMD::calcKineticPrim计算动能的Primitive估计值

$$T_{\mathrm{prim}} = \frac{1}{n} \sum_i^N \sum_j^n \frac{[p_i^{(j)}]^2}{2m_i} - \frac{1}{2n} \sum_i^N \sum_j^n m_i \omega_n^2 (x_i^{(j)} - x_i^{(j + 1)})^2 = T_p - T_q$$

(2) RPMD::calcKinetic计算动能的Virial估计值

$$T_v = \frac{1}{n^2} \sum_i^N \sum_j^n \frac{|\mathbf{p}_i^{(j)}|^2}{2m_i} - \frac{1}{2n} \sum_i^N \sum_j^n (\mathbf{x}_i^{(j)} - \overline{\mathbf{x}}_i) \cdot \mathbf{F}_i^{(j)} = T_p - T_q$$

(3) RPMD::calcKineticClassical计算经典动能估计值

$$T_c = \sum_i^N \sum_j^n \frac{[p_i^{(j)}]^2}{2m_i} + \frac{1}{2} \sum_i^N \sum_j^n m_i \omega_n^2 (x_i^{(j)} - x_i^{(j + 1)})^2$$

一般$T_c \neq nT_{\mathrm{prim}}$，NVE系综中的守恒量为$T_c$加上经典势能贡献

## 关于正则变换

In [1], the normal transformation is base on real fourier series, which is unusual in FFT libraries. Here we shall derive complex version of normal transformation.

The hamiltonian reads:

$$H = \sum_{j = 1}^n \frac{p_j^2}{2m} + \frac{1}{2} m \omega^2 (q_{j + 1} - q_j)^2$$

The complex fourier series reads:

$$\tilde{p}_k = \sum_{j = 1}^n p_j \exp(-i \frac{2 \pi}{n}jk)$$

$$p_j = \frac{1}{n} \sum_{k = 0}^{n - 1} \tilde{p}_k \exp(i \frac{2 \pi}{n}jk)$$

Substitude $p_j$ and $q_j$ using $\tilde{p}_k$ and $\tilde{q}_k$, we have

$$n^2 H = \sum_{j \alpha \beta} [\frac{1}{2m} \tilde{p}_\alpha \tilde{p}_\beta + \frac{1}{2} m \omega^2 \tilde{q}_\alpha \tilde{q}_\beta (e^{i \frac{2 \pi}{n} \alpha} - 1) (e^{i \frac{2 \pi}{n} \beta} - 1)] \exp(i \frac{2 \pi}{n} (\alpha + \beta))$$

$$n^2 H = \sum_{\alpha \beta} [\frac{1}{2m} \tilde{p}_\alpha \tilde{p}_\beta + \frac{1}{2} m \omega^2 \tilde{q}_\alpha \tilde{q}_\beta (e^{i \frac{2 \pi}{n} \alpha} - 1) (e^{i \frac{2 \pi}{n} \beta} - 1)] \sum_j \exp(i \frac{2 \pi}{n} (\alpha + \beta))$$

Making use of the fact that $\sum_j \exp(i \frac{2 \pi}{n} (\alpha + \beta)) = \delta_{\alpha + \beta, 0}$, we have

$$n^2 H = \sum_{\alpha} \frac{1}{2m} \tilde{p}_\alpha \tilde{p}_{-\alpha} + \frac{1}{2} m (4 \omega^2 \sin^2(\frac{\pi}{n} \alpha)) \tilde{q}_\alpha \tilde{q}_{-\alpha}$$

$$n^2 H = \sum_{\alpha} \frac{1}{2m} |\tilde{p}_\alpha|^2 + \frac{1}{2} m (4 \omega^2 \sin^2(\frac{\pi}{n} \alpha)) |\tilde{q}_\alpha|^2$$

$$n^2 H = \sum_{\alpha} \frac{1}{2m} |\tilde{p}_\alpha|^2 + \frac{1}{2} m \omega_\alpha^2 |\tilde{q}_\alpha|^2$$

where $\omega_\alpha = 2 \omega \sin(\frac{\pi}{n} \alpha))$

# 关于积分器

一般使用BAOAB和BCOCB积分器，演化算符

$$\hat S = \hat B \hat A \hat O \hat A \hat B$$

n个时间步后分子动力学的过程为$\hat S^n$，为了实现上的方便默认不初始化$t = 0$时的力，因此实际计算的是$\hat B^{-1} \hat S^n$。该区别对平衡态分子动力学影响有限。

## Reference

[1] J. Chem. Phys. 133, 124104 (2010); https://doi.org/10.1063/1.3489925
