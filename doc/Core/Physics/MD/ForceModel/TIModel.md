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

## 关于势能零点

记$U_0$为爱因斯坦晶体的势场，$U_1$为要研究体系的势场，考虑参数$\lambda$变化下混合系统的势场

$$U(\lambda) = \lambda U_1 + (1 - \lambda) U_0$$

其中$U_0$以完美晶格为势能零点，但$U_1$势能零点尚未指定，亥姆霍兹自由能$F$随参数变化

$$\langle \frac{\partial F}{\partial \lambda} \rangle = \langle \frac{\partial U(\lambda)}{\partial \lambda} \rangle = \langle U_1 - U_0 \rangle$$

$U_1$可相差一个任意常数，为消除该任意性，假定完美晶格为$U_1$的势能零点。TIModel::deltaPotentialV()的结果将是绝对的。

## Reference

[1] Comput. Mater. Sci. 112, 333-341 (2016); https://doi.org/10.1016/j.commatsci.2015.10.050
