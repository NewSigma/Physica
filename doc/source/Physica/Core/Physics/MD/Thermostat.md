# Defination of temperature

Derivation of $\tilde C_{\mathbf{vv}}(t)$

$$\tilde C_{\mathbf{rr}}(t) = \langle \frac{1}{N} \sum_i^N \mathbf{\overline r}_i(0) \cdot \mathbf{\overline r}_i(t) \rangle$$

$$\frac{\mathrm d}{\mathrm d t}C_{\mathbf{rr}}(t) = \langle \frac{1}{N} \sum_i^N \mathbf{\overline r}_i(0) \cdot \mathbf{\overline v}_i(t) \rangle = -\langle \frac{1}{N} \sum_i^N \mathbf{\overline v}_i(0) \cdot \mathbf{\overline r}_i(t) \rangle$$

$$\tilde C_{\mathbf{vv}}(t) = -\frac{\mathrm d^2}{\mathrm d t^2}C_{\mathbf{rr}}(t) = \langle \frac{1}{N} \sum_i^N \mathbf{\overline v}_i(0) \cdot \mathbf{\overline v}_i(t) \rangle = \langle \frac{1}{N} \sum_i^N \frac{\mathbf{\overline p}_i(0) \cdot \mathbf{\overline p}_i(t) \rangle}{m_i^2}$$

According to classical statistical physics

$$\tilde C_{\mathbf{vv}}(0) = \frac{3\sum_i^N m_i^{-1}}{N \beta}$$

Quantum estimator of temperature is defined by

$$T_{\mathrm{quantum}} = \frac{\sum_i^N \frac{|\mathbf{\overline p}_i(0)|^2}{m_i^2}}{3k \sum_i^N m_i^{-1}}$$

while classical estimator is computed using theorem of equipartition of energy

$$T_{\mathrm{classical}} = \frac{1}{3Nn^2k} \sum_i^N \sum_k^n \frac{|\mathbf p_i^{(k)}|^2}{m_i}$$

# Langevin Thermostat

$$\frac{\mathrm{d} p_j}{\mathrm{d}t} = -\gamma p_j(t) + \sqrt{\frac{2 m \gamma}{\beta_n}} \xi_j(t)$$

Using normal representation

$$p_j = \frac{1}{\sqrt n} \sum_{k = 0}^{n - 1} \tilde p_k e^{-i \frac{2\pi}{n}jk}$$

Langevin thermostat in normal representation

$$\frac{\mathrm{d} \tilde p_k}{\mathrm{d}t} = -\gamma \tilde p_k(t) + \sqrt{\frac{2 m \gamma}{\beta_n}} \tilde \xi_k(t)$$

where $\tilde \xi_k = \frac{1}{\sqrt n} \sum_{j = 1}^{n} \tilde \xi_k e^{-i \frac{2\pi}{n}jk}$

$$\mathbb{E}[\mathrm{Re}[\tilde \xi_k]] = 0 \qquad \mathbb{V}[\mathrm{Re}[\tilde \xi_k]] = \frac{1}{2}$$

where we have used

$$\sum_j \cos^2 \frac{2\pi}{n} jk = \frac{1}{2} + \frac{1}{2n} \sum_j \cos \frac{4\pi}{n} jk $$

$$\sum_j \cos \frac{4\pi}{n} jk = \sum_j \mathrm{Re}[\exp(i\frac{4\pi}{n} jk)] = \mathrm{Re}[\frac{1 - \exp(i 4\pi k)}{1 - \exp(i\frac{4\pi}{n} k)}] = 0$$

## Reference

[1] Miller TF, Manolopoulos DE. J. Chem. Phys. 2005. 122:184503
