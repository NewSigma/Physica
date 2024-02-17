# TIModel

## 爱因斯坦晶体的亥姆霍兹自由能

哈密顿量

$$H = \sum_i^N (\frac{\mathbf{p_i}^2}{2m_i} + \frac{1}{2} m_i \omega_i^2 \mathbf{r}_i^2)$$

爱因斯坦晶体不同粒子相互独立，因而配分函数是可分离的

$$Z = \prod_i^N Z_i$$

其中单粒子配分函数

$$Z_i = \frac{1}{h^3} \int \exp(\frac{-\mathbf{\beta p_i}^2}{2m_i}) \mathrm{d^3} p \int \exp(-\frac{\beta}{2} m_i \omega_i^2 \mathbf{r}_i^2) \mathrm{d^3} r$$

$$= \frac{1}{h^3} (\frac{2\pi m_i}{\beta})^{\frac 3 2} \frac{1}{\omega_i^3}(\frac{2\pi}{\beta m_i})^{\frac 3 2}$$

$$= (\frac{kT}{\hbar \omega_i})^3$$

亥姆霍兹自由能

$$F = -kT \ln Z = -kT \sum_i^N \ln Z_i = 3kT \sum_i^N \ln\frac{\hbar \omega_i}{kT}$$

若$\omega_i = \omega$，则$F = 3NkT \ln \frac{\hbar \omega}{kT}$，与文献[1]式(15)一致

## Reference

[1] Comput. Mater. Sci. 112, 333-341 (2016); https://doi.org/10.1016/j.commatsci.2015.10.050
