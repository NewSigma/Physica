# Fokker-Planck model

The implementation is devoted to solve 2D non-static Fokker-Planck function, consider dynamic equations:

$$\left \{ \begin{matrix} \dot{x} = \frac{p}{m} \\ \dot{p} = F(x, p) + \sqrt{2D(x, p)} \Gamma(t) \end{matrix} \right.$$

where $D(x, p)$ is diffusion coefficient, $\Gamma(t)$ is Gaussian white noise of unit magnitude.

The corresponding Fokker-Planck equation is

$$\frac{\partial W(x, p, t)}{\partial t} = [-\frac{p}{m} \frac{\partial}{\partial x} - \frac{\partial}{\partial p} F(x, p) - \frac{\partial}{\partial p} D \frac{\partial D}{\partial p} + \frac{\partial^2}{\partial p^2} D] W(x, p, t)$$

Note that $[\frac{\partial}{\partial p} D \frac{\partial D}{\partial p}]W = \frac{\partial}{\partial p} (WD \frac{\partial D}{\partial p})$

Let

$$x = \frac{\xi}{\sqrt{1 - \xi^2}} \qquad y = \frac{\eta}{\sqrt{1 - \eta^2}}$$

We have

$$\frac{\partial \tilde{W}(\xi, \eta, t)}{\partial t} = [-\frac{\eta}{m \sqrt{1 - \eta^2}} (1 - \xi^2)^{\frac{3}{2}} \frac{\partial}{\partial \xi} - (1 - \eta^2)^{\frac{3}{2}} \frac{\partial}{\partial \eta} F + (1 - \eta^2)^{\frac{3}{2}} \frac{\partial}{\partial \eta} ((1 - D) \frac{\partial D}{\partial p} + D (1 - \eta^2)^{\frac{3}{2}} \frac{\partial}{\partial \eta})] \tilde{W}(\xi, \eta, t)$$

To convert it into FEM format, we integrate and apply greens' formular, we obtain

$$\mathrm{RHS} = \sum_{kj} c_j \int [-\frac{\eta}{m \sqrt{1 - \eta^2}} (1 - \xi^2)^{\frac{3}{2}} u_i \frac{\partial u_j}{\partial \xi} + (F + (D - 1) \frac{\partial D}{\partial p}) \sqrt{1 - \eta^2} ((1 - \eta^2) \frac{\partial u_i} {\partial \eta} - 3 \eta u_i) u_j + 3D \eta (1 - \eta^2)^2 u_i \frac{\partial u_j}{\partial \eta} - D (1 - \eta^2)^3 \frac{\partial u_i} {\partial \eta} \frac{\partial u_j}{\partial \eta}] \mathrm{d} S_k$$

where $\tilde{W}(\xi, \eta, t) = \sum_i c_i(t) u_i(\xi, \eta)$
