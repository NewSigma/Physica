<!--
This file is part of Physica.

SPDX-FileCopyrightText: 2025-2026 Weibo He
SPDX-License-Identifier: GFDL-1.3-no-invariants-or-later
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
