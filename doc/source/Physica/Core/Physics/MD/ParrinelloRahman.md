# Algorithm for stress

Atom position in crystal: cartesian form $\mathbf{x}$ and direct form $\mathbf{s}$ is related by $\mathbf{x = As}$, where $\mathbf{A}$ is lattice matrix

Suppose lattice changed by $\delta \mathbf{A}$, displacement is $\delta \mathbf{x} = \delta \mathbf{A} \mathbf{A}^{-1} \mathbf{x}$

Cauchy strain is defined as $\mathbf{U}_{ij} = \frac{1}{2} (\frac{\partial \delta \mathbf{x}_i}{\partial \mathbf{x_j}} + \frac{\partial \delta \mathbf{x}_j}{\partial \mathbf{x_i}})$

(Green strain is defined as $\mathbf{U}_{ij} = \frac{1}{2} (\frac{\partial \delta \mathbf{x}_i}{\partial \mathbf{x_j}} + \frac{\partial \delta \mathbf{x}_j}{\partial \mathbf{x_i}} + \frac{\partial \delta \mathbf{x}^k}{\partial \mathbf{x_i}} \frac{\partial \delta \mathbf{x}_k}{\partial \mathbf{x_j}})$)

Therefore

$$\mathbf{U} = \frac{1}{2} (\delta \mathbf{A} \mathbf{A}^{-1} + \mathbf{A}^{-T} \delta \mathbf{A}^{T})$$

The Cauchy stress is defined as $\mathbf{S} = \frac{\partial E}{\partial \mathbf{U}}$, where $E = E(\mathbf{U(\delta A)})$ is internal energy.

$$\frac{\partial E}{\partial \delta \mathbf{A}^{m}_{\; n}} = \frac{\partial E}{\partial \mathbf{U}^{\alpha}_{\; \beta}} \frac{\partial \mathbf{U}^{\alpha}_{\; \beta}}{\partial \delta \mathbf{A}^i_{\; j}} = \mathbf{S}^{\beta}_{\;\alpha} \frac{\partial \mathbf{U}^{\alpha}_{\; \beta}}{\partial \delta \mathbf{A}^i_{\; j}}$$

We compute

$$\partial_{\delta \mathbf A} \mathbf{U}^{in}_{\quad jm} = \frac{\partial \mathbf U^i_{\; j}}{\partial \delta \mathbf{A}^{m}_{\; n}} = \frac{1}{2}[(\mathbf{A}^{-1})^n_{\; j} \delta^i_{\; m} + (\mathbf{A}^{-1})^n_{\; i} \delta^j_{\; m}]$$

Substitude into the equation above

$$(\partial_{\mathbf A} E)^{m}_{\; n} = (\mathbf{A}^{-1})^m_{\; \alpha} \mathbf{S}^\alpha_{\; n}$$

In matrix form

$$\mathbf{S} = \mathbf{A} (\partial_{\mathbf A} E)$$

Internal energy $E = T + V$

$$(\partial_{\mathbf A}T)^a_{\;b} = \frac{\partial T}{\partial \mathbf{A}^{b}_{\; a}} = \sum_i^N m_i \dot s^a_{\;(i)} v_{b(i)}$$

$$\mathbf{S}_\mathbf{T} = \mathbf{A}^m_{\;\alpha} (\partial_{\mathbf A}T)^\alpha_{\;n} = \sum_i^N m_i v^m_{\;(i)} v_{n(i)}$$

$$(\partial_{\mathbf A}V)^a_{\;b} = \frac{\partial V}{\partial \mathbf{A}^{b}_{\; a}} = \sum_i \frac{\partial V}{\partial \mathbf{r}^\alpha_{\;(i)}} \frac{\partial \mathbf{r}^\alpha_{\;(i)}}{\partial \mathbf{A}^{b}_{\; a}} = -\sum_i \mathbf{s}^a_{(i)} \mathbf{F}_{b(i)}$$

$$\mathbf{S}_\mathbf{V} = \mathbf{A}^m_{\;\alpha} (\partial_{\mathbf A}V)^\alpha_{\;n} = -\sum_i^N \mathbf{r}^m_{\;(i)} \mathbf{F}_{n(i)}$$

Cauchy stress is

$$\mathbf{S} = \sum_i^N m_i v^m_{\;(i)} v_{n(i)} -\sum_i^N \mathbf{r}^m_{\;(i)} \mathbf{F}_{n(i)}$$

由角动量守恒，$\mathbf S$是对称的

# Notes on P-R NPT ensemble

We shall follow notations in [1], lagrange equation reads

$$\frac{\partial L}{\partial s^a_{(i)}} = \frac{\mathrm d}{\mathrm d t} \frac{\partial L}{\partial \dot s^a_{(i)}}$$

$$\mathrm{LHS} = \mathbf{A}^\alpha_{\;a} F_{\alpha(i)}$$

$$\mathrm{RHS} = \frac{\mathrm d}{\mathrm d t} (\mathbf{A}^\alpha_{\;a} p_{\alpha(i)})$$

simplify

$$\mathbf{\dot p}^a_{(i)} = \mathbf{F}_{(i)} - \mathbf{A}^{-1} \mathbf{\dot A} \mathbf{p}$$

# Reference
[1] M. Parrinello and A. Rahman, J. Appl. Phys. 52, 7182 (1981); doi: 10.1063/1.328693
