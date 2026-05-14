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
# Vegas

## On Chi-Square Statistics

The chi-square statistic is used to evaluate the reliability of grid refinement$^{[1]}$:

$$\chi_0^2 = \sum_j \frac{(I_j - \overline{I})^2}{\sigma_j^2}$$

According to the central limit theorem, each term in the sum follows a normal distribution with standard deviation $(N_s)^{-1/2}$, where $N_s$ is the number of samples (NumSample). Therefore $N_s\chi_0^2$ follows a chi-square distribution with $N_r - 1$ degrees of freedom, where $N_r$ is the number of grid refinement iterations (NumRefine). Normalizing the above:

$$\chi^2 = \frac{N_s}{N_r - 1} \sum_j \frac{(I_j - \overline{I})^2}{\sigma_j^2}$$

Unlike the definition in [1], the above expression has an expected value of 1. Specifically, when $N_r = 1$, no grid refinement is involved, and we define $\chi^2 = 1$.

## Template Parameter

**TakeLn**: Averaging a series of samples spanning many orders of magnitude often leads to overflow or underflow. When `TakeLn = true`, the log-mean and log-variance are computed to alleviate this issue. The core idea is to use the identity $x = e^{\ln x} \quad (x > 0)$ to separate the largest magnitude component. Log-mean:

$$\ln\braket{x} = \ln\left( \frac{1}{n} \sum_i^n x_i \right) = \ln\left( \frac{1}{n} \sum_i^n e^{\ln x_i + \ln x_m - \ln x_m} \right) = \ln x_m + \ln\left( \frac{1}{n} \sum_i^n e^{\ln x_i - \ln x_m} \right)$$

$$= \ln x_m + \ln\braket{e^{\ln x_i - \ln x_m}}$$

Similarly, log-variance:

$$\ln \sigma^2(x) = \ln(\braket{x_i^2} - \braket{x_i}^2) = 2\ln x_m + \ln \sigma^2(e^{\ln x_i - \ln x_m})$$

where $x_m = \max(x_1, x_2, ..., x_n)$. Since $e^{\ln x_i - \ln x_m} \le 1$, summing these terms is generally well-behaved.

Weighted summation of results from different iteration counts yields statistics with lower uncertainty (Ref. [1] Eq. 30 ~ Eq. 32). The log-mode equivalents are:

$$\ln\overline{I} = \ln\sum_j e^{\ln I_j - \ln\sigma_j^2} + \ln\sigma_{\overline I}^2$$

$$\ln\sigma_{\overline I}^2 = -\ln \sum_j e^{-\ln\sigma_j^2}$$

$$\chi^2 = \sum_j\exp(2\ln\overline{I} + 2\ln|e^{\ln I_j - \ln\overline{I}} - 1| - \ln\sigma_j^2)$$

## Reference

[1] J. Comput. Phys. 439, 110386 (2021); https://doi.org/10.1016/j.jcp.2021.110386
