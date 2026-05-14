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
# TIModel

## Helmholtz Free Energy of the Einstein Crystal

Hamiltonian:

$$H = \sum_i^N (\frac{\mathbf{p_i}^2}{2m_i} + \frac{1}{2} m_i \omega_i^2 \mathbf{r}_i^2)$$

In an Einstein crystal, different particles are independent of each other, so the partition function is separable:

$$Z = \prod_i^N Z_i$$

The single-particle partition function is:

$$Z_i = \frac{1}{h^3} \int \exp(\frac{-\mathbf{\beta p_i}^2}{2m_i}) \mathrm{d^3} p \int \exp(-\frac{\beta}{2} m_i \omega_i^2 \mathbf{r}_i^2) \mathrm{d^3} r$$

$$= \frac{1}{h^3} (\frac{2\pi m_i}{\beta})^{\frac 3 2} \frac{1}{\omega_i^3}(\frac{2\pi}{\beta m_i})^{\frac 3 2}$$

$$= (\frac{kT}{\hbar \omega_i})^3$$

Helmholtz free energy:

$$F = -kT \ln Z = -kT \sum_i^N \ln Z_i = 3kT \sum_i^N \ln\frac{\hbar \omega_i}{kT}$$

If $\omega_i = \omega$, then $F = 3NkT \ln \frac{\hbar \omega}{kT}$, consistent with Eq. (15) of Ref. [1].

## On the Zero Point of Potential Energy

Let $U_0$ be the potential field of the Einstein crystal, and $U_1$ be the potential field of the system under study. Consider the mixed system potential under variation of the parameter $\lambda$:

$$U(\lambda) = \lambda U_1 + (1 - \lambda) U_0$$

where $U_0$ takes the perfect lattice as its zero of potential energy, but the zero of $U_1$ has not yet been specified. The Helmholtz free energy $F$ varies with the parameter as:

$$\langle \frac{\partial F}{\partial \lambda} \rangle = \langle \frac{\partial U(\lambda)}{\partial \lambda} \rangle = \langle U_1 - U_0 \rangle$$

$U_1$ can differ by an arbitrary constant. To eliminate this arbitrariness, we assume the perfect lattice is the zero of potential energy for $U_1$. The result of `TIModel::deltaPotentialV()` will then be absolute.

## Reference

[1] Comput. Mater. Sci. 112, 333-341 (2016); https://doi.org/10.1016/j.commatsci.2015.10.050
