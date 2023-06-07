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

## Reference

[1] Miller TF, Manolopoulos DE. J. Chem. Phys. 2005. 122:184503
