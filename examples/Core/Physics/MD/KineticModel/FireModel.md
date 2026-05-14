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
# FIRE Structure Optimization

## Introduction

Atomic structure optimization is an important step in first-principles simulations and molecular dynamics simulations. The FIRE structure optimization algorithm has significant advantages over traditional conjugate gradient methods, and its performance can even rival more complex quasi-Newton methods$^{[1]}$. This document demonstrates how to use the FIRE algorithm for fixed-volume and fixed-pressure structure optimization. We also discovered that the fixed-timestep FIRE algorithm outperforms the adaptive-timestep version, providing evidence to guide algorithm selection in practical applications.

## Fixed-Volume Structure Optimization (FireModel)

![](./FireModel.png)

**Figure 1** After adding a Gaussian perturbation to the optimized structure and optimizing with FIRE, the force magnitude as a function of iteration step converges to machine epsilon for double precision after 1500 steps.

### Fixed-Timestep FIRE

The FIRE algorithm$^{[1,2]}$ uses an initial timestep $\Delta t$ with an adaptive timestep $\tau_n$. The semi-implicit Euler algorithm of FIRE$^{[2]}$ can be written as:

$$\mathbf{v}'_{n + 1} = \mathbf{v}_n + \frac{\mathbf F}{m} \tau_n$$

$$\mathbf{v}_{n + 1} = (1 - \alpha) \mathbf{v}'_{n + 1} + \alpha \frac{|\mathbf v'_{n + 1}|}{|\mathbf F|} \mathbf{F}$$

$$\mathbf{x}_{n + 1} = \mathbf{x}_n + \mathbf{v}_{n + 1} \tau_n$$

where the second step ensures, according to the vector collinearity theorem, that the velocity magnitude decreases and the system energy is lowered.

The fixed-timestep version shifts the adaptive part onto $\lambda_n$:

$$\mathbf{v}'_{n + 1} = \mathbf{v}_n + \lambda_n \frac{\mathbf F}{m} \Delta t$$

$$\mathbf{v}_{n + 1} = (1 - \alpha) \mathbf{v}'_{n + 1} + \alpha \frac{|\mathbf v'_{n + 1}|}{|\mathbf F|} \mathbf{F}$$

$$\mathbf{x}_{n + 1} = \mathbf{x}_n + \mathbf{v}_{n + 1} \Delta t$$

where $\lambda_n = \tau_n / \Delta t$, meaning the fixed-timestep version can be implemented using adaptive-timestep dynamics steps. Since the mixing step only requires the unit vector of the force, this modification has no effect on the mixing step. The fixed-timestep version effectively adaptively adjusts the force magnitude, and experiments show this modification significantly improves optimization efficiency.

Attempting to reproduce the optimization of molten SiO2 at 6000K from Ref. [2]:

![](./SiO2_1.png)

**Figure 1** Adaptive-timestep FIRE algorithm optimizing molten SiO2 structure: force magnitude vs. iteration count

![](./SiO2_2.png)

**Figure 2** Fixed-timestep FIRE algorithm optimizing molten SiO2 structure: force magnitude vs. iteration count

From these two figures, the fixed-timestep version clearly outperforms the adaptive-timestep version. However, the adaptive-timestep version still surpasses the literature results for reasons unknown.

## Fixed-Pressure Structure Optimization (CFireModel)

The original FIRE algorithm$^{[1,2]}$ does not consider lattice degrees of freedom. The FIRE implementations in mainstream MD and first-principles software such as LAMMPS$^{[4]}$ and QE$^{[5]}$ do not support variable-cell optimization. Some literature combines position and lattice degrees of freedom but does not disclose the details$^{[6]}$.

Starting from the anisotropic SCR dynamics equation, which is equivalent to the Parrinello-Rahman equation$^{[3]}$ in the overdamped limit, we convert the first-order equation to second-order:

$$\left\{ \begin{matrix} \mathrm{d} \mathbf{h} = \frac{\mathbf{\Pi}}{W} \mathrm{d}t \\ \mathrm{d} \mathbf{\Pi} = -\mathbf{Dh} \mathrm{d}t \end{matrix} \right.$$

where $\mathbf{\Pi}$ is the generalized momentum of the lattice degrees of freedom, and $\mathbf{G = -Dh}$ is the generalized force. For convergence determination, the FIRE algorithm uses the 2-norm of the force; here, lattice convergence similarly uses the 2-norm of the generalized force as the criterion. Integration similarly uses the fixed-timestep semi-implicit Euler method, and the mixing step is:

$$\mathbf{\Pi} = (1 - \alpha) \mathbf{\Pi}' + \alpha \frac{|\mathbf{\Pi}|}{|\mathbf{G}|}\mathbf{G}$$

![](./CFireModel1.png)

**Figure 1** Force and generalized force vs. iteration count. After 1500 iterations, both force and generalized force are reduced to machine precision.

![](./CFireModel2.png)

**Figure 2** Pressure vs. iteration count. After 1500 iterations, the pressure stabilizes at $10^{-18}$, which is smaller than double-precision machine epsilon, and thus is not a convergence criterion at machine precision.

**Assertion**: The default value of the lattice mass can be set to the number of free components in the lattice matrix.

## 与其他软件对比

![](./QE/OptQE.png)

**Figure 1** Force vs. iteration count. The fixed-timestep FIRE algorithm converges faster than QE's adaptive-timestep FIRE. FIRE1 (mass not normalized) takes about 250 steps to reach the target precision, while FIRE2 (mass normalized) takes about 125 steps.

![](./QE/OptQE1.png)

**Figure 2** Force vs. iteration count. With a fixed 200 steps, Fire achieves higher precision than Damping.

## Reference

[1] Phys. Rev. Lett. 97, 170201 (2006); https://doi.org/10.1103/PhysRevLett.97.170201
[2] Comput. Mater. Sci. 175, 109584 (2020); https://doi.org/10.1016/j.commatsci.2020.109584  
[3] https://doi.org/10.48550/arXiv.2111.06402  
[4] LAMMPS Doc; https://docs.lammps.org/min_style.html  
[5] QE Doc; https://www.quantum-espresso.org/Doc/INPUT_PW.html  
[6] Phys. Rev. B 93, 020508(R); https://doi.org/10.1103/PhysRevB.93.020508  
