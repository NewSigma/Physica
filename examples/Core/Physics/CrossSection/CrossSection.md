# Some notes on implementation

The radial Schrdinger equation reads:

$$\frac{d^2 u(r)}{d r^2} = [\frac{2m}{\hbar^2} (V(r) - E) + \frac{l(l + 1)}{r^2}]u(r)$$

In implementation we let $h = \frac{r}{\rho}$, so the equation becomes

$$\frac{d^2 U(h)}{d h^2} = [\frac{2m \rho^2}{\hbar^2} (V(r) - E) + \frac{l(l + 1)}{h^2}]U(h)$$

Here $\frac{2m \rho^2}{\hbar^2}$ is the source of $6.12 meV^{-1} \rho^{-2}$ provided in [1].

Reference:
[1] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013:21
