<!--
This file is part of Physica.

SPDX-FileCopyrightText: 2025-2026 Weibo He
SPDX-License-Identifier: GFDL-1.3-no-invariants-or-later
-->

# DiffBerry - AutoDiff approach to Berry curvature

This example uses automatic differentiation to compute the Berry curvature and band Chern number. The advantage of automatic differentiation is that it requires only a single eigenstate, providing an elegant approach to circumvent the well-known gauge ambiguity inherent in finite difference methods.

Consider a two-band Chern insulator on a two-dimensional lattice $^{[1]}$:

$$H(\mathbf k) = \vec{d}(\mathbf{k}) \cdot \vec{\sigma}$$

where $\vec{d}(\mathbf{k}) = (\sin(k_x), 3\sin(k_y), 1 - \cos(k_x) - \cos(k_y))$, $\vec{\sigma} = (\sigma_x, \sigma_y, \sigma_z)$ is Pauli's matrix. The Berry curvature reads

$$\Omega(\mathbf{k}) = -2\; \text{Im}[\braket{\frac{\partial \psi}{\partial k_x}|\frac{\partial \psi}{\partial k_y}}]$$

Here $\ket{\psi}$ is a nondegenerate eigenstate. The Chern numbers of the ground state and the excited state are $-1$ and $1$, respectively.

## Reference

[1] Phys. Rev. B 91, 165140; <https://doi.org/10.1103/PhysRevB.91.165140>
