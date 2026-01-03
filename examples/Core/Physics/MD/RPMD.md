<!--
Copyright 2024-2026 Weibo He.

This file is part of Physica.

Permission is granted to copy, distribute and/or modify this document
under the terms of the GNU Free Documentation License, Version 1.3
or any later version published by the Free Software Foundation;
with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.

You should have received a copy of the GNU Free Documentation License
along with Physica.  If not, see <https://www.gnu.org/licenses/>.
-->
# RPMD

## Introduction

在分子动力学模拟中，我们通常采用波恩奥本海默近似。即考虑到原子核质量远大于电子质量，忽略电子和声子的相互作用。原子核坐标作为体系哈密顿量的参数，其运动规律由牛顿第二定律描述。但是在含有轻元素的体系中，该近似不能很好的解释实验现象。比如室温下水的RDF:

![](./Analyser/RDF.png)

**图1** 298K下水的氢原子间RDF，PIMD结果与文献[1]吻合良好。由于核量子效应PIMD的第一个峰显著低于MD的第一个峰。

在经典分子动力学(MD)模拟中，由于未考虑核量子效应(NQE)，MD的第一个峰显著高于考虑了核量子效应的路径积分分子动力学(PIMD)的第一个峰。对比文献[1]发现，PIMD的结果更接近实验。可以得出结论：在常温常压下水的核量子效应是显著的。理论分析指出，随着温度降低或压强增加，一些较重的元素也会表现出显著的核量子效应。模拟核量子效应的数值方法主要有路径积分蒙特卡罗(PIMC)，路径积分分子动力学和量子热浴方法(QTB)。其中，PIMC不能计算系统的动力学性质。而QTB相比PIMD忽略了不同弹簧环间声子气体的散射效应，所以一般应以PIMD为基准$[2]$。因此本文只讨论PIMD方法。

## PIMD

PIMD可以视作MD的推广，当副本(Replica)数$n = 1$时(详见下文)，PIMD将退化为MD。注意在讨论PIMD时“经典”可以有电子和原子核两方面的含义:

|          | MD | AIMD | PIMD | AI-PIMD |
| -------- | -- | ---- | ---- | ------- |
| 电子 | 经典 | 量子 | 经典 | 量子 |
| 原子核 | 经典 | 经典 | 量子 | 量子 |

MD相当于在玻尔兹曼分布中进行采样，对N粒子系统，正则系综配分函数为

$$Q_c = \frac{1}{(2\pi\hbar)^N} \int e^{-\beta H(\mathbf{p, q})} \text{d}\mathbf{p}\text{d}\mathbf{q}$$

类似地, PIMD相当于在量子玻尔兹曼分布中进行采样。将积分改为求迹，将哈密顿量改为哈密顿算符可得(不严格地):

$$Q = \text{tr}[e^{-\beta \hat{H}}]$$

迹一般是难以计算的，我们利用路径积分近似计算配分函数。将路径分成$n$段:

$$Q = \text{tr}[e^{-\beta \hat{H}}] = \text{tr}[(e^{-\frac{\beta}{n} \hat{H}})^n]$$

利用Trotter分解，由于动能和势能算符一般不对易，简单拆开将导致小的误差

$$e^{-\frac{\beta}{n} \hat{H}} = e^{-\frac{\beta}{n} \hat{T}} e^{-\frac{\beta}{n} \hat{V}} + o(\frac{1}{n^2})$$

显然$n$越大，该近似越精确。丢弃二阶小量，配分函数近似为

$$Q \approx \text{tr}[\Lambda_0 e^{-\frac{\beta}{n} \hat{T}} e^{-\frac{\beta}{n} \hat{V}} \lambda_1 e^{-\frac{\beta}{n} \hat{T}} e^{-\frac{\beta}{n} \hat{V}} \Lambda_2 ... \Lambda_{n - 1} e^{-\frac{\beta}{n} \hat{T}} e^{-\frac{\beta}{n} \hat{V}} \Lambda_n]$$

在每个$\Lambda_i$处插入一组正交完备的基矢

$$\text{id} = \int \text{d} q_i \text{d} p_i \ket{q_i} \braket{q_i|p_i} \bra{p_i}$$

我们熟知从路径积分的定义出发进行计算是繁琐的，感兴趣的读者可参考相关文献$^{[3]}$。我们只给出化简结果并着重讨论它的物理意义:

$$Q \approx \frac{1}{(2\pi\hbar)^{Nn}} \int e^{-\beta_n H_n} \text{d}\mathbf{q} \text{d}\mathbf{p}$$

其中

$$H_n = \sum^N_{i = 1} \sum^n_{j = 1} [\frac{[p_i^{(j)}]^2}{2m} + \frac{1}{2}m\omega_n^2 (q_i^{(j)} - q_i^{(j + 1)})^2 + V(q_i^{(j)})]; \quad \beta_n = \beta/n \quad \omega_n = (\beta_n \hbar)^{-1}$$

当$n = 1$时，$H_1$是经典系统的哈密顿量，对应着PIMD退化为MD。当$n \neq 1$, 上式相当于将粒子数为$N$, 温度为$\beta$的量子系统映射为粒子数为$nN$, 温度为$\beta_n$的经典系统。为了理解$H_n$的物理意义，我们注意到第二项为简谐势贡献，其弹性常数$k \propto \omega_n^2 = (\beta_n \hbar)^{-2}$包含量子力学的特征常数$\hbar$。因此第二项为纯粹的量子效应。若不考虑量子贡献，配分函数可以分离变量

$$Q = Q_c^n$$

上式表明, 如不考虑耦合项，$H_n$将分成$n$个完全独立的经典系统。量子效应将这$n$个经典系统用等效的弹簧耦合。

为了理解这个图像，我们设想存在$n$个平行宇宙，每个宇宙都存在一个装有经典理想气体的容器。每个容器内理想气体具有相同的粒子数$N$、体积$V$和温度$T$。我们对每个粒子给定序号$i \quad (1 \le i \le N)$。我们对所有粒子进行如下操作：按照 $1, 2, 3, ..., (n - 1), n, 1$的顺序将 $n$个宇宙的 $i$粒子用弹性系数为 $k$的弹簧首尾相连，形成一个弹簧环。由于处在不同的宇宙当中，除了弹簧外，不同宇宙的粒子不存在静电、引力、交换等相互作用。这里的弹簧正是核量子效应。现在我们将 $n$个平行宇宙“拍”到一个宇宙当中，我们就得到一个相对符合实际的宇宙。即考虑到量子统计效应后，任何物体都是由$n$个经典对象组成的一个弹簧环，物体不再是质点而是存在一个确定的密度分布。这个图像为我们理解量子力学和量子统计提供了一个有趣的角度。

一般来说$\hbar$是一个相对很小的量，这导致以$\hbar^2$为分母的弹性常数$k$很大，弹簧近乎刚性。因此被弹簧连接的粒子将局域在一个很小的空间范围内。由于 $k \propto T^2$，温度越高弹性系数越大，核量子效应越不显著。在高温情况下，可以近似认为弹簧几乎不会伸长，弹簧环退化为质点，量子将过渡到经典情况。

在理解了PIMD的物理图像后，我们重新回顾上述理论框架。在量子力学的路径积分表述中，从初态到末态的所有可能路径都使用一个作用量的指数因子$\exp(\frac{i}{\hbar} S)$进行加权求和。但是在上面的推导中，我们完全没有用到系统的作用量 $S$。因此PIMD中的路径积分并非“物理”的路径积分，而单纯只是一种计算配分函数的数学技巧。它只能给出近似的量子动力学。在前文讲述的平行宇宙图像中，我们对所有粒子进行了序号标记。这意味着我们先验地假设了粒子是可分辨的，即忽略了粒子间的交换相互作用。这个近似在温度不十分低的凝聚态体系中可以适用，但是在温度达到十几K甚至更低的时候就必须考虑粒子间的交换相互作用。近年来费米子和玻色子的PIMD理论也得到一些发展，感兴趣的读者可以关注。

## 结语

最后引用文献[2]中的一段话:

    If one performs only classical simulations,
    one will never know whether quantum effects are important.
    One must have the ability to include quantum effects into a simulation,
    even if only approximately,
    to know when they are important and when they are not.

先验地判断一个体系中核量子效应是否重要是困难的。在下定结论前可以使用PIMD等方法进行简单测试，有可能会发现新奇的现象。

## Reference

[1] J. Chem. Phys. 131, 024501 (2009); https://doi.org/10.1063/1.3167790
[2] Rev. Mod. Phys. 89, 035003 (2017); https://doi.org/10.1103/RevModPhys.89.035003
[3] A. Altland, B. Simons Condensed Matter Field Theory 2nd edition[M]. Cambridge University Press, 2010:97-100
