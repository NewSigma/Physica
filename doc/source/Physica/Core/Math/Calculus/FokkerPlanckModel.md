# Fokker-Planck model

The implementation is devoted to solve 2D non-static Fokker-Planck function, consider dynamic equations:

$$\left \{ \begin{matrix} \dot{x} = \frac{p}{m} \\ \dot{p} = F(x, p) + \sqrt{2D(x, p)} \Gamma(t) \end{matrix} \right.$$

where $D(x, p)$ is diffusion coefficient, $\Gamma(t)$ is Gaussian white noise of unit magnitude.

The corresponding Fokker-Planck equation is

$$\frac{\partial W(x, p, t)}{\partial t} = [-\frac{p}{m} \frac{\partial}{\partial x} - \frac{\partial}{\partial p} F(x, p) - \frac{\partial}{\partial p} D \frac{\partial D}{\partial p} + \frac{\partial^2}{\partial p^2} D] W(x, p, t)$$

Note that $[\frac{\partial}{\partial p} D \frac{\partial D}{\partial p}]W = \frac{\partial}{\partial p} (WD \frac{\partial D}{\partial p})$

To convert it into FEM format, we integrate and apply greens' formular, we obtain

$$\mathrm{RHS} = \sum_{kj} c_j \int \{ -\frac{p}{m} u_i \frac{\partial u_j}{\partial x} + [F + (D - 1) \frac{\partial D}{\partial p}] u_i  \frac{\partial u_i}{\partial p} - D \frac{\partial u_i}{\partial p} \frac{\partial u_j}{\partial p} \} \mathrm{d} S_k$$

where $\tilde{W}(\xi, \eta, t) = \sum_i c_i(t) u_i(\xi, \eta)$
