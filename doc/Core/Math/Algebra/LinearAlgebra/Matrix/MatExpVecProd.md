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
# MatExpVecProd

Calculates matrix exponential-vector product$^{[1]}$:

$$\mathbf{y} = e^{\mathbf{A}} \mathbf{x}$$

## Template param

**NoFactor**: Section 3.1 of [1] mentioned ways to approximate $e^{\mathbf{A}}$:

$$e^\mu[T_m(\frac{1}{n} (A -\mu I))]^n$$

If the constant $e^\mu$ is not essential, we may ignore it, which will improve numerical stability and performance.

## Reference

[1] SIAM J. Sci. Comput. 33(2), 488–511 (2011); <https://doi.org/10.1137/100788860>
