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
# TPQ

This example demonstrates how to prepare TPQ states and compute observables using TPQ states.

## Introduction

Thermo Pure Quantum (TPQ) states$^{[1]}$ are used to approximately compute various mechanical and thermodynamic quantities of finite-temperature many-body systems. As the number of lattice sites increases, the computed results converge exponentially to the exact solution in a probabilistic sense. This convergence is quite fast: even with 4 lattice sites, TPQ results are already very close to exact diagonalization results. The key improvement of TPQ over exact diagonalization is that it reduces the required $O(N)$ eigenvectors for finite-temperature calculations to a single random vector, dramatically reducing memory requirements.

## Numerical Results

![](./TPQ.png)

**Figure 1** Comparison of *Physica* and *HPhi* 3.5.2$^{[2]}$ for computing the relationship between inverse temperature $\beta$ and particle density in the Hubbard model, using both full diagonalization and TPQ methods, with $U/t = 8$, lattice sites $L = 4$. The statistical uncertainty is smaller than the line width. Since HPhi cannot adaptively adjust the expansion order, its results accumulate numerical error more quickly.

    # HPhi Input (FullDiag)
    model = "Fermion HubbardGC"
    method = "Full Diag"
    lattice = "chain"
    L = 4
    U = 8
    t = 1

    # HPhi Input (TPQ)
    model = "Fermion HubbardGC"
    method = "cTPQ"
    lattice = "chain"
    L = 4
    U = 8
    t = 1
    lanczos_max = 161
    LargeValue = 10
    NumAve = 8192
    # Increase expansion order to reduce deviation between TPQ and Full-ED
    # ExpandCoef = 10

Benchmarking on an Intel(R) Xeon(R) Platinum 8358 + 256G platform, *Physica* is 1300x faster than HPhi 3.5.2 for a 4x4 Hubbard model.

## Reference

[1] Phys. Rev. Lett. 111, 010401 (2013); <https://doi.org/10.1103/PhysRevLett.108.240401>
[2] HPhi; <https://github.com/issp-center-dev/HPhi.git>
