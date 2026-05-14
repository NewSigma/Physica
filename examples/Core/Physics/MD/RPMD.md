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
# RPMD

## Introduction

In molecular dynamics simulations, we typically employ the Born-Oppenheimer approximation, which assumes that because the mass of atomic nuclei is much larger than that of electrons, the interaction between electrons and phonons can be neglected. The nuclear coordinates serve as parameters in the system Hamiltonian, and their motion is described by Newton's second law. However, in systems containing light elements, this approximation cannot adequately explain experimental observations. For example, the RDF of water at room temperature:

![](./Analyser/RDF.png)

**Fig. 1** RDF between hydrogen atoms in water at 298K. The PIMD results agree well with Ref. [1]. Due to nuclear quantum effects, the first peak of PIMD is significantly lower than that of MD.

In classical molecular dynamics (MD) simulations, because nuclear quantum effects (NQE) are not considered, the first peak of MD is significantly higher than that of path integral molecular dynamics (PIMD), which does account for NQE. Comparison with Ref. [1] shows that PIMD results are closer to experiment. This leads to the conclusion that nuclear quantum effects in water are significant at ambient temperature and pressure. Theoretical analysis indicates that as temperature decreases or pressure increases, even heavier elements exhibit significant nuclear quantum effects. The main numerical methods for simulating NQE include path integral Monte Carlo (PIMC), path integral molecular dynamics (PIMD), and the quantum thermal bath (QTB) method. PIMC cannot compute dynamical properties of the system. QTB neglects the phonon gas scattering effects between different ring polymers compared to PIMD, so PIMD should generally serve as the benchmark$[2]$. Therefore, this document only discusses the PIMD method.

## PIMD

PIMD can be seen as a generalization of MD. When the number of replicas $n = 1$ (detailed below), PIMD reduces to MD. Note that when discussing PIMD, "classical" can refer to both electrons and nuclei:

|          | MD | AIMD | PIMD | AI-PIMD |
| -------- | -- | ---- | ---- | ------- |
| Electrons| Classical | Quantum | Classical | Quantum |
| Nuclei   | Classical | Classical | Quantum | Quantum |

MD corresponds to sampling from the Boltzmann distribution. For an N-particle system, the canonical ensemble partition function is:

$$Q_c = \frac{1}{(2\pi\hbar)^N} \int e^{-\beta H(\mathbf{p, q})} \text{d}\mathbf{p}\text{d}\mathbf{q}$$

Similarly, PIMD corresponds to sampling from the quantum Boltzmann distribution. Replacing the integral with a trace and the Hamiltonian with the Hamiltonian operator (informally):

$$Q = \text{tr}[e^{-\beta \hat{H}}]$$

The trace is generally difficult to compute. We use the path integral approximation to compute the partition function by splitting the path into $n$ segments:

$$Q = \text{tr}[e^{-\beta \hat{H}}] = \text{tr}[(e^{-\frac{\beta}{n} \hat{H}})^n]$$

Using the Trotter decomposition, since kinetic and potential energy operators generally do not commute, simply splitting them introduces a small error:

$$e^{-\frac{\beta}{n} \hat{H}} = e^{-\frac{\beta}{n} \hat{T}} e^{-\frac{\beta}{n} \hat{V}} + o(\frac{1}{n^2})$$

Clearly, the larger $n$ is, the more accurate this approximation. Dropping the second-order term, the partition function is approximately:

$$Q \approx \text{tr}[\Lambda_0 e^{-\frac{\beta}{n} \hat{T}} e^{-\frac{\beta}{n} \hat{V}} \lambda_1 e^{-\frac{\beta}{n} \hat{T}} e^{-\frac{\beta}{n} \hat{V}} \Lambda_2 ... \Lambda_{n - 1} e^{-\frac{\beta}{n} \hat{T}} e^{-\frac{\beta}{n} \hat{V}} \Lambda_n]$$

Inserting a complete set of basis vectors at each $\Lambda_i$:

$$\text{id} = \int \text{d} q_i \text{d} p_i \ket{q_i} \braket{q_i|p_i} \bra{p_i}$$

Carrying out the full derivation from the path integral definition is tedious; interested readers may refer to the relevant literature$^{[3]}$. We only give the simplified result and focus on its physical interpretation:

$$Q \approx \frac{1}{(2\pi\hbar)^{Nn}} \int e^{-\beta_n H_n} \text{d}\mathbf{q} \text{d}\mathbf{p}$$

where

$$H_n = \sum^N_{i = 1} \sum^n_{j = 1} [\frac{[p_i^{(j)}]^2}{2m} + \frac{1}{2}m\omega_n^2 (q_i^{(j)} - q_i^{(j + 1)})^2 + V(q_i^{(j)})]; \quad \beta_n = \beta/n \quad \omega_n = (\beta_n \hbar)^{-1}$$

When $n = 1$, $H_1$ is the Hamiltonian of the classical system, corresponding to PIMD reducing to MD. When $n \neq 1$, the above formula maps a quantum system with $N$ particles at temperature $\beta$ to a classical system with $nN$ particles at temperature $\beta_n$. To understand the physical meaning of $H_n$, note that the second term is a harmonic potential contribution whose spring constant $k \propto \omega_n^2 = (\beta_n \hbar)^{-2}$ contains the characteristic quantum constant $\hbar$. Therefore, the second term is a purely quantum effect. Without quantum contributions, the partition function separates variables:

$$Q = Q_c^n$$

This shows that without the coupling term, $H_n$ splits into $n$ completely independent classical systems. Quantum effects couple these $n$ classical systems through effective springs.

To visualize this, imagine there exist $n$ parallel universes, each containing a container filled with a classical ideal gas. Each container has the same number of particles $N$, volume $V$, and temperature $T$. We assign an index $i \quad (1 \le i \le N)$ to each particle. We then perform the following operation on all particles: connect the $i$-th particles from $n$ universes end-to-end in the order $1, 2, 3, ..., (n - 1), n, 1$ using a spring of stiffness $k$, forming a ring. Since they exist in different universes, particles from different universes have no electrostatic, gravitational, exchange, or other interactions except through the springs. These springs represent nuclear quantum effects. Now, if we "collapse" the $n$ parallel universes into one, we obtain a relatively realistic universe. That is, after accounting for quantum statistical effects, any object consists of a spring ring made of $n$ classical copies. Objects are no longer point masses but have a well-defined density distribution. This picture provides an interesting perspective for understanding quantum mechanics and quantum statistics.

In general, $\hbar$ is relatively small, which makes the spring constant $k$, with $\hbar^2$ in the denominator, very large -- the springs are nearly rigid. Consequently, the particles connected by springs are localized within a very small spatial region. Since $k \propto T^2$, higher temperatures lead to larger spring constants and less significant nuclear quantum effects. At high temperatures, the springs can be approximately considered inextensible, the spring ring degenerates into a point particle, and the quantum case transitions to the classical one.

Having understood the physical picture of PIMD, we revisit the theoretical framework above. In the path integral formulation of quantum mechanics, all possible paths from the initial to the final state are weighted by an exponential factor of the action $\exp(\frac{i}{\hbar} S)$. However, in the derivation above, we never used the action $S$ of the system. Thus, the path integral in PIMD is not a "physical" path integral; it is purely a mathematical technique for computing the partition function. It can only provide approximate quantum dynamics. In the parallel universe picture described earlier, we labeled all particles with indices. This means we a priori assumed particles are distinguishable, i.e., we neglected exchange interactions between particles. This approximation is valid in condensed matter systems at not-too-low temperatures, but at temperatures of tens of K or lower, exchange interactions must be considered. In recent years, PIMD theory for fermions and bosons has also seen some development; interested readers may follow up.

## Conclusion

Finally, quoting from Ref. [2]:

    If one performs only classical simulations,
    one will never know whether quantum effects are important.
    One must have the ability to include quantum effects into a simulation,
    even if only approximately,
    to know when they are important and when they are not.

It is difficult to determine a priori whether nuclear quantum effects are important in a system. Before drawing conclusions, one can perform simple tests using methods like PIMD, and may discover novel phenomena.

## Reference

[1] J. Chem. Phys. 131, 024501 (2009); https://doi.org/10.1063/1.3167790
[2] Rev. Mod. Phys. 89, 035003 (2017); https://doi.org/10.1103/RevModPhys.89.035003
[3] A. Altland, B. Simons Condensed Matter Field Theory 2nd edition[M]. Cambridge University Press, 2010:97-100
