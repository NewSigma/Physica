# Rotation of spherical hamonic function

## Real spherical hamonic function

The spherical hamonic function in the code is defined as following:

$$Y_{lm}(\theta, \phi) = \sqrt{\frac{2l + 1}{2 \pi \delta_m} \frac{(l - |m|)!}{(l + |m)!|}} P^{|m|}_l(\cos{\theta}) \Phi(\phi)$$

where $\delta_m$ and $\Phi(\phi)$ is defined as following:

$$
\delta_m =
\left \{
    \begin{matrix}
        2 \quad (m = 0) \\
        1 \quad (m \neq 0)
    \end{matrix}
\right.
$$

$$
\Phi(\phi) = 
\left \{
    \begin{matrix}
        \cos(m \phi) \quad (m \geq 0) \\
        \sin(m \phi) \quad (m < 0)
    \end{matrix}
\right.
$$

The rotation matrix $R$ acts like this

$$\langle \hat Y| = \langle Y|R$$

Relation between the defination and defination of complex function is

$$Y^C_{lm} = \sqrt{\frac{\delta_m}{2}} Y^R_{l,|m|} + i\sqrt{\frac{\delta_m}{2}} Y^R_{l,-|m|}$$

where $Y^C_{lm}$ is complex spherical hamonic function and $Y^R_{lm}$ is real spherical hamonic function.
