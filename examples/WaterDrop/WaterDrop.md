# Surface of the water drop

## Basic theory

Suppose we have a water drop on a infinite smooth plain.

The water drop has axial symmetry so the surface of the water drop can be described using function $z(r)$. Build up cylinder coordinate system whose z-axis is along the axis of symmetry of the drop. The infinite plain is at $z = 0$.

We have the initial condition

$$z'(0) = 0 \quad z(r_0) = 0$$

where $r_0$ is radius of the water drop.

Energy of the system consists of gravitational potential energy $E_g$ and surface energy $E_s$, that is

$$E = E_g + E_s$$

$$E_g = \int g z(r) dm = \pi \rho g \int^{r_0}_0 r [z(r)]^2 dr$$

$$E_s = \sigma(S_1 + S_2) = \sigma \pi {r_0}^2 + 2 \pi \sigma \int_0^{r_0} \sqrt{1 + [z'(r)]^2} r dr$$

where $S_1$ is the square of contact surface between the drop and plain, $S_2$ is the square of surface of the drop.

Because the conservation of mass, we have

$$\int dm = 2 \pi \rho \int_0^{r_0} r z(r) dr = M$$

The functional of energy is

$$J = E[z(r)] + \lambda [2 \pi \rho \int_0^{r_0} r z(r) dr - M]$$

$$\delta J = \int_0^{r_0} [\rho g r z(r) - \sigma \frac{r z''(r) + z'(r) + [z'(r)]^3}{(1 + [z'(r)]^2)^{\frac{3}{2}}} + \lambda \rho r] \delta z dr = 0$$

where we have made use of $\delta z|_{r = 0} = \delta z|_{r = r_0} = 0$.

We can optain

$$\sigma \frac{r z''(r) + z'(r) + [z'(r)]^3}{(1 + [z'(r)]^2)^{\frac{3}{2}}} = \rho g r z(r) + \lambda \rho r$$

The initial value of the ODE is

$$
\left \{
    \begin{matrix}
        z(r_0) = 0 \\
        z'(r_0) = -tan(\theta)
    \end{matrix}
\right.
$$

where $\theta$ is contact angle.

The solution is

$$z = z(r, r_0, \lambda)$$

But we have two constraint equations:

$$M = \int dm$$

$$z'(0) = 0$$

## Numerical method

Let $f(r) = z(-r)$, substitude into the ODE, we obtain

$$\sigma \frac{r f''(r) + f'(r) + [f'(r)]^3}{(1 + [f'(r)]^2)^{\frac{3}{2}}} = \rho g r f(r) + \lambda \rho r$$

The initial value of is

$$f(-r_0) = 0 \quad f'(-r_0) = tan(\theta)$$

Here we can use forward Runge-Kutta method. Let $p = \frac{df(r)}{dr}$, we obtain

$$
\left \{
    \begin{matrix} \frac{df}{dr} = p \\
    \frac{dp}{dr} = \frac{(1 + p^2)^{\frac{3}{2}} (\rho g f(r) + \lambda \rho)}{\sigma} - \frac{p(1 + p^2)}{r}
    \end{matrix}
\right.
$$

$$f(-r_0) = 0 \quad p(-r_0) = tan(\theta)$$
