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

General molecular dynamics simulations consist of four core components:

    KineticModel: equations of motion, integration algorithm, boundary conditions, configurational constraints
    ForceModel: force field
    Thermostat: temperature control algorithm, NVT, NPT
    Barostat: pressure control algorithm, NPT only

These four components are organically integrated by the class RPMD. Taking a 1D system as an example, the RPMD Hamiltonian is:

$$H_n = \sum_i^N \sum_j^n [\frac{[p_i^{(j)}]^2}{2m_i} + \frac{1}{2} m_i \omega_n^2 (q_i^{(j)} - q_i^{(j + 1)})^2] + \sum_j^n V(\mathbf{q}^{(j)})$$

where $n$ is the number of replicas (NumReplica). When $n = 1$, the Hamiltonian reduces to that of classical molecular dynamics (MD). The class RPMD aims to implement MD and RPMD within a unified framework.

The partition function under the canonical ensemble with fixed boundary conditions is:

$$Q = \frac{1}{(2\pi \hbar)^{Nn}} \int e^{-\beta_n H_n(\mathbf{p}, \mathbf{q})} \mathrm{d} \mathbf{p} \mathrm{d} \mathbf{q}$$

The partition function under the canonical ensemble with periodic boundary conditions is:

$$Q = \frac{1}{(2\pi \hbar)^{Nn}} \int e^{-\beta_n H_n(\mathbf{p}, \mathbf{q})} \delta(\sum_{ij} p_i^{(j)}) \mathrm{d} \mathbf{p} \mathrm{d} \mathbf{q}$$

## Kinetic Energy Calculation

(1) `RPMD::calcKineticPrim` computes the primitive estimator of kinetic energy:

$$T_{\mathrm{prim}} = \frac{1}{n} \sum_i^N \sum_j^n \frac{[p_i^{(j)}]^2}{2m_i} - \frac{1}{2n} \sum_i^N \sum_j^n m_i \omega_n^2 (x_i^{(j)} - x_i^{(j + 1)})^2 = T_p - T_q$$

(2) `RPMD::calcKinetic` computes the virial estimator of kinetic energy:

$$T_v = \frac{1}{n^2} \sum_i^N \sum_j^n \frac{|\mathbf{p}_i^{(j)}|^2}{2m_i} - \frac{1}{2n} \sum_i^N \sum_j^n (\mathbf{x}_i^{(j)} - \overline{\mathbf{x}}_i) \cdot \mathbf{F}_i^{(j)} = T_p - T_q$$

(3) `RPMD::calcKineticClassical` computes the classical estimator of kinetic energy:

$$T_c = \sum_i^N \sum_j^n \frac{[p_i^{(j)}]^2}{2m_i} + \frac{1}{2} \sum_i^N \sum_j^n m_i \omega_n^2 (x_i^{(j)} - x_i^{(j + 1)})^2$$

In general $T_c \neq nT_{\mathrm{prim}}$. The conserved quantity in the NVE ensemble is $T_c$ plus the classical potential energy contribution.

## On Normal Transformation

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

## On Integrators

The BAOAB and BCOCB integrators are commonly used. The evolution operator is:

$$\hat S = \hat B \hat A \hat O \hat A \hat B$$

After $n$ time steps, the molecular dynamics process is $\hat S^n$. For implementation convenience, the force at $t = 0$ is not initialized by default, so the actual computation is $\hat B^{-1} \hat S^n$. This difference has a limited effect on equilibrium molecular dynamics.

## Reference

[1] J. Chem. Phys. 133, 124104 (2010); https://doi.org/10.1063/1.3489925
