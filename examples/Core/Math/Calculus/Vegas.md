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
# Vegas - 自适应蒙特卡洛积分

Vegas算法自适应地调整积分网格以降低采样误差, 考虑8维积分

$$I = \int_{-a}^a \text{d}\mathbf{x} \exp(-\frac{1}{2} \mathbf{x}^T \mathbf{x})$$

定义网格损失函数为所有维度中$d_i$标准差的最大值

$$L = \text{max}_i \; \sigma(d_i) = \text{max}_i \sqrt{\frac{\braket{d_i^2} - \braket{d_i}^2}{\braket{d_i}^2}}$$

其中$d_i$由[1]中(17)式定义.

![](./Vegas1.png)

**图1** 网格损失随迭代次数变化. compressRate控制网格调整的速度, 压缩率越大网格收敛越快, 过大的压缩率将导致不稳定$^{[1]}$

![](./Vegas2.png)

**图2** 更密的网格通常能够达到更好的网格损失

![](./Vegas3.png)

**图3** 样本数不足将导致网格改善过程中的不稳定

可以使用卡方统计评估网格改善的可靠性, $\chi^2$的理论值为$1$, 显著大于1可能意味着遗失了部分被积函数的细节.

## Reference

[1] J. Comput. Phys. 439, 110386 (2021); https://doi.org/10.1016/j.jcp.2021.110386
