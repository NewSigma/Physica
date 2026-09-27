<!--
Copyright 2025-2026 Weibo He.

This file is part of Physica.

Permission is granted to copy, distribute and/or modify this document
under the terms of the GNU Free Documentation License, Version 1.3
or any later version published by the Free Software Foundation;
with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.

You should have received a copy of the GNU Free Documentation License
along with Physica.  If not, see <https://www.gnu.org/licenses/>.
-->
# Notes on Implementation of DQMC

## Rank-1 Update

We apply the $\Delta$ matrix from the right, so the update formulas differ slightly from [1]. Ignoring spin indices, Eq. (7.43) and Eq. (7.45) should be modified respectively to:

$$R = 1 + (1 - G_{ii})\Delta_{ii}$$

$$G \to G - \frac{1}{R}(I - G)\Delta G$$

## Relative sign

In quantum Monte Carlo, the expectation value of an observable $\hat A$ is calculated using absolute-value reweighting:

$$\braket{\hat A} = \frac{\sum_i A_i |w_i| s_i}{\sum_i |w_i| s_i},$$

where $\sum_i$ runs over all possible configurations, and $w_i = |w_i| s_i$ is the weight of configuration $i$. In DQMC, $s_i$ may be obtained from a recursive relation $s_{i + 1} = f(s_i)$, which is cheaper than computing it from scratch, while the approach still requires knowing $s_0$. Note that the numerator and denominator may differ by a coefficient without changing $\braket{\hat A}$. Let **relative sign** $r_i = \frac{s_i}{s_0} $. Obviously, $r_0 = 1$, and the recursive relation still holds. We obtain

$$\braket{\hat A} = \frac{\sum_i A_i |w_i| s_0 r_i}{\sum_i |w_i| s_0 r_i} = \frac{\sum_i A_i |w_i| r_i}{\sum_i |w_i| r_i}.$$

Thus we avoid computing $s_0$ if we do not actually need it.

## Reference

[1] Gubernatis J, Kawashima N, Werner P. Quantum Monte Carlo Methods: Algorithms for Lattice Models. Cambridge University Press; 2016:194  
