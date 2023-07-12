# About translational invariance

两个问题:
1. 由于有限的计算精度，合力在三个方向(x, y, z)上的分量均不等于0 ($\Gamma$点上声学支频率不为0)
2. 在简谐近似下，系统势能

$$V \approx V_0 + \frac 1 2 \sum_{im\alpha, jn\beta} \Phi_{im\alpha, jn\beta} u_{im\alpha} u_{jn\beta}$$

力常数有对称性

$$\Phi_{im\alpha, jn\beta} = \frac{\partial^2 V}{\partial u_{im\alpha} \partial u_{jn\beta}} = \frac{\partial^2 V}{\partial u_{jn\beta} \partial u_{im\alpha}} = \Phi_{jn\beta, im\alpha}$$

简谐振动下，力常数可由下式计算

$$\Phi_{im\alpha, jn\beta} = -\frac{F_{im\alpha, jn\beta}}{u_{jn\beta}}$$

由于非谐效应，上式是近似的，得到的力常数不满足对称性，力常数矩阵不对称，特征向量不是正交的。仿照[1]使用迭代方法解决

## Interpolation of force constant matrix

For simplicity, we write force constants in k-space

$$\tilde \Phi_{mn}(k_n) = \frac{1}{N_c} \sum_{i = 1}^{N_c} \Phi_{mn}(i) e^{-j k_n a_i}$$

where $N_c$ is the number of cell, $j = \sqrt{-1}$ is imaginary unit, $a_i = ia$ is the position of i th unit cell in r-space, $k_n = \frac{n}{N_c} \frac{2\pi}{a} \quad (n = 0, 1, ..., N_c - 1)$ is position of k points, $k_n a_i = \frac{2\pi in}{N_c}$

使用傅里叶基内插倒空间力常数矩阵，需要相距更远的两原子间的力常数，由于力随距离衰减截断力常数可作为一种选择，但绝不是唯一的

$$\Phi_{mn}(i) = \left\{ \begin{matrix} \neq 0 \quad (1 \leq i \leq N_c) \\ = 0 \quad (N_c < i \leq N_c') \end{matrix} \right.$$

内插的倒空间力常数矩阵

$$\tilde \Phi_{mn}(k_n) = \frac{1}{N_c'} \sum_{i = 1}^{N_c'} \Phi_{mn}(i) e^{-j k_n a_i}$$

$k_n = \frac{n}{N_c'} \frac{2\pi}{a} \quad (n = 0, 1, ..., N_c' - 1)$, 该矩阵有$N_c' > N_c$个k点

## Reference

[1] Dario Alfè PHON: A program to calculate phonons using the small displacement method [J]. Computer Physics Communications, 2009, 180(12), 2622-2633
