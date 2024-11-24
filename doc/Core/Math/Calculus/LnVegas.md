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
# LnVegas

## Introduction

对一系列大数量级样本进行平均通常导致上溢, 可以通过计算对数均值和对数方差规避该问题. 核心思想是利用恒等式$x = e^{\ln x} \quad (x > 0)$, 将数量级最大的部分分离出去

对数均值

$$\ln\braket{x} = \ln(\frac{1}{n} \sum_i^n x_i) = \ln(\frac{1}{n} \sum_i^n e^{\ln x_i + \ln x_m - \ln x_m}) = \ln x_m + \ln(\frac{1}{n} \sum_i^n e^{\ln x_i - \ln x_m}) = \ln x_m + \ln\braket{e^{\ln x_i - \ln x_m}}$$

同理, 对数方差

$$\ln \sigma^2(x) = \ln(\braket{x_i^2} - \braket{x_i}^2) = 2\ln x_m + \ln \sigma^2(e^{\ln x_i - \ln x_m})$$

其中$x_m = \max(x_1, x_2, ..., x_n)$. 由于$e^{\ln x_i - \ln x_m} \le 1$, 对其求和的性质是良好的.

该技术对下溢同样适用

## Vegas变式

使用上述方法确保Vegas积分中的稳定性, 同样可以对不同迭代次数的结果进行加权求和以得到不确定度更低的统计量([1]中Eq. 30 ~ Eq. 32):

$$\ln\overline{I} = \ln\sum_j e^{\ln I_j - \ln\sigma_j^2} + \ln\sigma_{\overline I}^2$$

$$\ln\sigma_{\overline I}^2 = -\ln \sum_j e^{-\ln\sigma_j^2}$$

$$\chi^2 = \sum_j\exp(2\ln\overline{I} + 2\ln|e^{\ln I_j - \ln\overline{I}} - 1| - \ln\sigma_j^2)$$

## Reference

[1] J. Comput. Phys. 439, 110386 (2021); https://doi.org/10.1016/j.jcp.2021.110386
