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
# Some notes on implementation

The radial Schrdinger equation reads:

$$\frac{d^2 u(r)}{d r^2} = [\frac{2m}{\hbar^2} (V(r) - E) + \frac{l(l + 1)}{r^2}]u(r)$$

In implementation we let $h = \frac{r}{\rho}$, so the equation becomes

$$\frac{d^2 U(h)}{d h^2} = [\frac{2m \rho^2}{\hbar^2} (V(r) - E) + \frac{l(l + 1)}{h^2}]U(h)$$

Here $\frac{2m \rho^2}{\hbar^2}$ is the source of $6.12 meV^{-1} \rho^{-2}$ provided in [1].

## Reference

[1] J. H. Thijssen. Computational Physics[M]. London: Cambridge University Press, 2013:21
