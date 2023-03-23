# Notes on P-R NPT ensemble

We shall follow notations in [1].

Equations of motion is

$$\ddot s_i = h^{-1} \sum_{j \neq i} \frac{F_{ij}}{m_i} - G^{-1} \dot G \dot s_i$$

or

$$\dot p_i = F_i - h^{-T} \dot G h^{-1} p_i$$

$$\dot p_i = F_i - \frac{1}{\Omega} (\dot h \sigma^T + (\dot h \sigma^T)^T) p_i$$

$$\dot p_i = F_i - \frac{1}{2 \pi} (\dot h B^T + (\dot h B^T)^T) p_i$$

Descrete

$$p_i := ((I - \frac{1}{2 \pi} (\dot h B^T + (\dot h B^T)^T)) p_i + F_i) \Delta t$$

where $p_i$ is momentum of i th particle and $b$ is matrix of reciprocal lattice.

$$\dot G = \frac{\partial G_{ij}}{\partial t}$$

# Reference
[1] M. Parrinello and A. Rahman, J. Appl. Phys. 52, 7182 (1981); doi: 10.1063/1.328693
