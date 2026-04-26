<!--
Copyright 2024-2026 Weibo He.

This file is part of Physica.

Permission is granted to copy, distribute and/or modify this document
under the terms of the GNU Free Documentation License, Version 1.3
or any later version published by the Free Software Foundation;
with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.

You should have received a copy of the GNU Free Documentation License
along with Physica.  If not, see <https://www.gnu.org/licenses/>.
-->
# JacobiDavidson

## Correction equation

The Jacobi–Davidson method gradually improves the approximate eigenvector by solving the correction equation

$$(\mathbf{I - UU}^H)(\mathbf A - \tilde \lambda \mathbf I) (\mathbf{I - UU}^H) \mathbf t = \mathbf{Bt = -r}$$

**Note**:

1. In practice, $\mathbf r \perp \mathbf U$ should be maintained; otherwise, the matrix $\mathbf B$ will be ill-conditioned and lead to numerical instability. That is, the modulus of the solution vector becomes too large, causing the orthogonalization to fail, which intuitively results in a pseudo-degeneracy phenomenon.

2. Since the new eigenvector will be sought in the subspace $\{ \mathbf{U, t} \}$, the sign of $\mathbf r$ in the above expression is not important.

## Implementation of refinedSearch()

The refined Ritz searching $^{[1]}$ requires finding $\mathbf x$ such that

$$\hat{\mathbf x} = \argmin_{|\mathbf x| = 1} |(\mathbf{AU} - \lambda \mathbf U)\mathbf x|$$

Let $\mathbf{B = AU- \lambda U}$, accordingly

$$|(\mathbf{AU - \lambda U})\mathbf x| = \sqrt{\mathbf x^H \mathbf B^H \mathbf B \mathbf x}$$

Note that $\mathbf B^H \mathbf B$ is positive definite. The minimum value of the above expression is exactly the minimum eigenvalue of $\mathbf B^H \mathbf B$, achieved when $\mathbf x$ is the eigenvector corresponding to that minimum eigenvalue. The optimization problem is thus reduced to an eigenvalue problem for a small matrix.

## Reference

[1] GAMM-Mitteilungen 29(2), 368-382 (2006); https://doi.org/10.1002/gamm.201490038
