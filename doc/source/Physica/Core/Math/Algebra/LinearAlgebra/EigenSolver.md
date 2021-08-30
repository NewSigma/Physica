# How to solve upper quasi-triangular linear system

If we want to compute the eigenvectors from real schur decomposition, we have to solve a upper quasi-triangular equation.

The algorithm solving this equation is clarified as following:

## $\lambda$ is real
The equation to solve[1] is:

$$(T_{11} - \lambda I)w = -u$$

The equation can be written as

$$\left[ \begin{matrix} T & u \end{matrix} \right] \left[ \begin{matrix} w \\ 1 \end{matrix} \right] = 0$$

The last element of solution vector is known as 1. Suppose we have got the value from $w_{i + 1}$ to $w_{n}$ we want to calculate $w_i$
(1) If ${T_{ii-1}} = 0$, so we have

$$(w_i, w_{i + 1}, ..., w_n, 1) \cdot (T_{ii}, ..., T_{in}, u_i) = 0$$

$$\implies w_i = \frac{-(w_{i + 1}, ..., w_n, 1) \cdot (T_{ii + 1}, ..., T_{in}, u_i)}{T_{ii}}$$

Note that $u_{i + 1} ... u_n$ is not used to compute $w_i$, so it is efficient to store $w_i$ to $u_i$

(2) If ${T_{ii-1}} < 0$, we have

$$\left[ \begin{matrix} T_{i - 1\ i - 1} & ... & T_{i - 1\ n} & u_{i - 1} \\ T_{i\ i - 1} & ... & T_{i\ n} & u_{i} \end{matrix} \right] \left[ \begin{matrix} w_{i - 1} \\ w_i \\ ... \\ w_n \\ 1 \end{matrix} \right] = 0$$

Let

$$a = -(w_{i + 1}, ..., w_n, 1) \cdot (T_{i - 1, i + 1}, ..., T_{i - 1, n}, u_{i - 1})$$

$$b = -(w_{i + 1}, ..., w_n, 1) \cdot (T_{i, i + 1}, ..., T_{i, n}, u_{i})$$

The equation becomes to

$$\left[ \begin{matrix} T_{i - 1\ i - 1} & T_{i - 1\ i} \\ T_{i\ i - 1} & T_{i\ i} \end{matrix} \right] \left[ \begin{matrix} w_{i - 1} \\ w_i \end{matrix} \right] = \left[ \begin{matrix} a \\ b \end{matrix} \right]$$

Let's use clammer's rule to solve it. The determinate of the coefficient matrix $C$ is

$$det(C) = \left| \begin{matrix} T_{i - 1\ i - 1} & T_{i - 1\ i} \\ T_{i\ i - 1} & T_{i\ i} \end{matrix} \right| = (Re(\lambda_C)) - \lambda)^2 + Im(\lambda_C)^2$$

where $\lambda_C$ is one of eigenvalues of $C$.

## $\lambda$ is complex
The equation to solve is

$$ \left[ \begin{matrix} T_{11} & u_1 & u_2 & T_{14} \\ 0 & a & b & v^T_1 \\ 0 & c & d & v^T_2 \\ 0 & 0 & 0 & T_{41} \end{matrix} \right] \left[ \begin{matrix} w \\ w' \\ w'' \\ 0 \end{matrix} \right] = \lambda \left[ \begin{matrix} w \\ w' \\ w'' \\ 0 \end{matrix} \right] $$

where value of $w'$ and $w''$ can be fetched by computing the eigenvectors of $\left[ \begin{matrix} a & b \\ c & d \end{matrix} \right]$, we have

$$w' = \frac{\frac{a - d}{2} \pm i \sqrt{|\frac{(a - d)^2}{4} + bc|}}{c}$$

$$w'' = 1$$

Now the equation becomes to

$$(T_{11} - \lambda I)w + w' u_1 + u_2 = 0$$

Complex eigenvalues corresponds to complex eigenvectors, so we have to divide the equation into two parts for the convenience of implementation. Note that $w''$ is a real number, so it is necessary to initialize vector $r$ on the right hand side of vector $c$ and the matrix becomes to a upper triangular matrix.

Let

$$w = r + ic$$

$$A = T_{11} - Re(\lambda)I$$

$$B = -Im(\lambda)I$$

The equation becomes to

$$\left \{ \begin{matrix} Ar - Bc + Re(w') u_1 + u_2 = 0 \\ Ac + Br + Im(w') u_1 = 0 \end{matrix} \right.$$

Suppose we have got values of $r_{i + 1} ... r_{n}$ and $c_{i + 1} ... c_{n}$ we want to calculate $w_i$

(1) If $A_{i i - 1} = 0$, we have

$$(A_{ii}, A_{i, i + 1}, ..., A_{in}, u_1, u_2) \cdot (r_i, r_{i + 1}, ..., r_n, Re(w'), 1) - B_{ii} c_i = 0$$

$$(A_{ii}, A_{i, i + 1}, ..., A_{in}, u_1, u_2) \cdot (c_i, c_{i + 1}, ..., c_n, Im(w'), 0) + B_{ii} r_i = 0$$

Let

$$\xi = -(A_{i, i + 1}, ..., A_{i, n}, u_1, u_2) \cdot (r_{i + 1}, ..., r_n, Re(w'), 1)$$

$$\eta = -(A_{i, i + 1}, ..., A_{i, n}, u_1, u_2) \cdot (c_{i + 1}, ..., c_n, Im(w'), 0)$$

We can get the solution:

$$ \left \{ \begin{matrix} r_i = \frac{A_{ii} \xi + B_{ii} \eta}{A_{ii}^2 + B_{ii}^2} \\ c_i = \frac{A_{ii} \eta - B_{ii} \xi}{A_{ii}^2 + B_{ii}^2} \end{matrix} \right. $$

Analogitically, it is efficient to save $r_i$ and $c_i$ to $u_1$ and $u_2$.

(2) If $A_{i i - 1} < 0$, we have

$$(A_{ii - 1}, A_{ii}, ..., A_{in}, u_1, u_2) \cdot (r_{i - 1}, r_{i}, ..., r_n, Re(w'), 1) - B_{ii} c_i = 0$$

$$(A_{ii - 1}, A_{i, i}, ..., A_{in}, u_1, u_2) \cdot (c_{i - 1}, c_{i}, ..., c_n, Im(w'), 0) + B_{ii} r_i = 0$$

$$(A_{i - 1 i - 1}, A_{i - 1 i}, ..., A_{i - 1 n}, u_1, u_2) \cdot (r_{i - 1}, r_{i}, ..., r_n, Re(w'), 1) - B_{i - 1 i - 1} c_{i - 1} = 0$$

$$(A_{i - 1 i - 1}, A_{i - 1, i}, ..., A_{i - 1 n}, u_1, u_2) \cdot (c_{i - 1}, c_{i}, ..., c_n, Im(w'), 0) + B_{i - 1 i - 1} r_{i - 1} = 0$$

Let

$$\xi_1 = -(A_{i, i + 1}, ..., A_{i, n}, u_1, u_2) \cdot (r_{i + 1}, ..., r_n, Re(w'), 1)$$

$$\xi_2 = -(A_{i, i + 1}, ..., A_{i, n}, u_1, u_2) \cdot (c_{i + 1}, ..., c_n, Im(w'), 0)$$

$$\xi_3 = -(A_{i - 1, i + 1}, ..., A_{i - 1, n}, u_1, u_2) \cdot (r_{i + 1}, ..., r_n, Re(w'), 1)$$

$$\xi_4 = -(A_{i - 1, i + 1}, ..., A_{i - 1, n}, u_1, u_2) \cdot (c_{i + 1}, ..., c_n, Im(w'), 0)$$

The equation can be written as

$$\left[ \begin{matrix} A_{i i - 1} & 0 & A_{ii} & -B_{ii} \\ 0 & A_{i i - 1} & B_{ii} & A_{ii} \\ A_{i - 1 i - 1} & -B_{i - 1 i - 1} & A_{i - 1 i} & 0 \\ B_{i - 1 i - 1} & A_{i - 1 i - 1} & 0 & A_{i - 1 i} \end{matrix} \right] \left[ \begin{matrix} r_{i - 1} \\ c_{i - 1} \\ r_i \\ c_i \end{matrix} \right] = \left[ \begin{matrix} \xi_1 \\ \xi_2 \\ \xi_3 \\ \xi_4 \end{matrix} \right]$$

Solution of this equation is complex, refer to [2] for a neat numerical algorithm.

## References:
[1] Golub, GeneH. Matrix computations = 矩阵计算 / 4th edition[M]. 人民邮电出版社, 2014.
[2] Eigen https://eigen.tuxfamily.org/
