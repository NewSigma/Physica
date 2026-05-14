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
# JacobiDavidson

## Introduction

The large-scale sparse eigenvalue problem (SEP) has extremely important and widespread applications in science and engineering. Numerical algorithms for solving SEP mainly include the power method$^{[1]}$, Lanczos$^{[1]}$, Arnoldi$^{[2]}$, Davidson$^{[3]}$, Jacobi-Davidson$^{[4]}$, FEAST$^{[5]}$, and LOBPCG$^{[6]}$. The power method can only solve for extreme eigenvalues, limiting its applicability. The Lanczos algorithm can be viewed as a special case of the Arnoldi algorithm for Hermitian matrices. The Arnoldi algorithm is used by mainstream scientific computing software such as Matlab$^{[7]}$ and Mathematica$^{[8]}$. Lanczos and Arnoldi algorithms, which are improvements on the power method, encounter difficulties with interior eigenvalue problems. The Davidson algorithm is fast for diagonally dominant matrices, but struggles to converge for non-diagonally dominant matrices. The Jacobi-Davidson algorithm introduces Jacobi orthogonal complement correction into the Davidson algorithm, maintaining Davidson's efficiency while solving its limitation to diagonally dominant matrices$^{[4]}$. The FEAST algorithm solves for all eigenvalues within a given interval, but is slower than Jacobi-Davidson. FEAST is implemented in Mathematica. LOBPCG is not suitable for non-symmetric matrices. The authors of both LOBPCG and JDCG (the symmetric variant of Jacobi-Davidson) claim their algorithm outperforms the other$^{[6, 10]}$, which is understandably confusing.

The Jacobi-Davidson algorithm has the following advantages:

- Fast convergence: Jacobi-Davidson converges faster than Davidson, Lanczos, and FEAST$^{[4, 9]}$
- Wide applicability: Jacobi-Davidson handles both extreme and interior eigenvalues. Unlike Davidson and LOBPCG, it imposes no requirements on the matrix type
- Strong modularity: Jacobi-Davidson decomposes into several independent components that can be freely combined to suit the problem

## Main

Consider the eigenvalue problem

$$\mathbf{Ax} = \lambda \mathbf{x}$$

where the eigenvector $\mathbf{x}$ and eigenvalue $\lambda$ are to be solved. We can randomly generate a vector $\mathbf{x_0}$ and substitute it into the above equation; the equality will not hold, giving a residual

$$\mathbf{r}_0(\lambda) = \mathbf{Ax}_0 - \lambda \mathbf{x}_0 \neq \mathbf{0}$$

At this point $\lambda$ is still unknown. We can do better than random guessing by using the Rayleigh quotient, which is the optimal estimate of the eigenvalue in the least-squares sense:

$$\theta_0 = \argmin_\lambda |\mathbf{r}_0(\lambda)| = \frac{\mathbf{x}_0^T \mathbf{Ax}}{\mathbf{x}_0^T \mathbf{x}}$$

Clearly $\mathbf{r}_0(\theta_0)$ is orthogonal to $\mathbf{x_0}$. The random initial guess $\mathbf{x_0}$ is not satisfactory, so we use Jacobi orthogonal complement correction to iteratively improve the solution accuracy.

Let the exact solution be expressed as corrections: $\mathbf{x} = \mathbf{x}_k + \delta \mathbf{x}_k$, $\lambda = \theta_k + \delta \theta_k$, where $\mathbf{x}_k, \theta_k$ are the eigenvector and eigenvalue at the $k$-th iteration. Substituting into the eigenvalue equation:

$$\mathbf{A} (\mathbf{x}_k + \delta \mathbf{x}_k) = (\theta_k + \delta \theta_k) (\mathbf{x}_k + \delta \mathbf{x}_k)$$

Rearranging and neglecting the second-order term $\delta \theta_k \cdot \delta \mathbf{x}_k$:

$$(\mathbf{A} - \theta I) \delta \mathbf{x}_k - \mathbf{x}_k \delta \theta_k = -\mathbf{r}_k(\theta_k)$$

To solve this equation, we multiply both sides by a projection matrix $I - \mathbf{x}_k \mathbf{x}_k^T$; the second term vanishes after projection:

$$(I - \mathbf{x}_k \mathbf{x}_k^T) (\mathbf{A} - \theta I) \delta \mathbf{x}_k = -\mathbf{r}_k(\theta_k)$$

Clearly, corrections along the $\mathbf{x}_k$ direction do not help improve the approximate solution. The $\delta \mathbf{x}_k$ of interest must be orthogonal to $\mathbf{x}_k$. Restricting the solution space to the orthogonal complement of $\mathbf{x}_k$:

$$(I - \mathbf{x}_k \mathbf{x}_k^T) (\mathbf{A} - \theta I) (I - \mathbf{x}_k \mathbf{x}_k^T) \delta \mathbf{x}_k = -\mathbf{r}_k(\theta_k)$$

This is the Jacobi-Davidson correction equation. The correction proceeds iteratively until convergence:

$$\mathbf{x}_{k + 1} = \mathbf{x}_k + \delta \mathbf{x}_k$$

$$\theta_{k + 1} = \argmin_\lambda |\mathbf{r}_{k + 1}(\lambda)|$$

In physics, we are typically interested in the ground state and the lowest few excited states, i.e., the smallest eigenvalues of a Hermitian matrix $\mathbf{A}$. Y. Notay proposed a variant of the Jacobi-Davidson algorithm for this case, using the conjugate gradient (CG) method to solve the Jacobi-Davidson correction equation. If the matrix is not Hermitian or if non-extreme eigenvalues are of interest, the linear solver module can be switched to more general solvers like MINRES or GMRES. This modular replaceability is one of the advantages of the Jacobi-Davidson algorithm.

## Reference

[1] Golub, GeneH. Matrix computations = 矩阵计算 / 4th edition[M]. 人民邮电出版社, 2014.284-285  
[2] Quart. Appl. Math. 9, 17–29 (1951); https://doi.org/10.1090/QAM/42792  
[3] J. Comput. Phys. 17, 87–94 (1975); https://doi.org/10.1016/0021-9991(75)90065-0  
[4] SIAM Review 42(2), 267–293 (2000); https://doi.org/10.1137/S0036144599363084  
[5] Phys. Rev. B 79, 115112 (2009); https://doi.org/10.1103/PhysRevB.79.115112  
[6] SIAM J. Sci. Comput. 23(2) 517–541 (2001); https://doi.org/10.1137/S1064827500366124  
[7] https://ww2.mathworks.cn/help/matlab/ref/eigs.html  
[8] https://reference.wolfram.com/language/ref/Eigenvalues.html  
[9] J. Comput. Electron. 14, 593–603 (2015); https://doi.org/10.1007/s10825-015-0695-z  
[10] Numerical Linear Algebra with Applications 9(1), 21-44 (2001); https://doi.org/10.1002/nla.246  
