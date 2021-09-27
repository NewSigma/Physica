# Notes on HFSolver

## Nomalization of wave function

In a general HF procedure, we require the wave function is normalized:

$$c^T S c = 1$$

where $c$ is the wave function vector, $S$ is the overlap matrix.

The passage aim to clarify the normalization is unnecessary in the code.

### Prove

When we solve the generalized eigenvalue problem:

$$F c = \epsilon S c$$

It is a usual method to turn it into a normal eigenvalue problem:

$$F' c' = \epsilon c'$$

There is relation between $c$ and $c'$:

$$c' = L c$$

where $L$ is a orthogonal matrix.

The resulting $c'$ is normalized using EigenSolver:

$$c'^T c' = 1$$

So we have

$$c^T S c = c'^T L S L^T c' = c'^T I c' = c'^T c' = 1$$

## Properties of density matrix

The density matrix $D$ is defined as

$$D = \sum^{N}_{i = 1} \alpha_i c_i c_i^T$$

where $N$ is the number of orbitals.
$\alpha_i = 1\ or\ 2$ is the number of electrons on orbital $i$, which satisfies $\sum^N_{i = 1} \alpha_i = N_e$, where $N_e$ is the total number of electrons.

The density matrix has the following properties:

(1) $\qquad Tr(SD) = N_e$
(2) $\qquad DSD = 3D - 2 \hat{D} \ge D$

where $\hat{D}$ is density matrix of electrons with up spin, comparation between matrices is element-wise.

### Prove

(1) $\qquad Tr(SD) = N_e$

$$Tr(SD) = Tr(L^{-1}SDL) = Tr(L^{-1}LL^TDL) = Tr(L^TDL)$$

where

$$Tr(L^T D L) = Tr(L^T \sum^N_{i = 1} \alpha_i c_i c_i^T L) = \sum^N_{i = 1} \alpha_i Tr(c_i' c'_{i} \ ^T) = \sum^N_{i = 1} \alpha_i |c_i|^2 = \sum^N_{i = 1} \alpha_i = N_e$$

(2) $\qquad DSD = 3D - 2 \hat{D} \ge D$

$$DSD = \sum_{i, j} \alpha_i \alpha_j c_i c_i^T L L^T c_j c_j^T$$

$$= \sum_{i, j} \alpha_i \alpha_j c_i c_i'^T c_j' c_j^T$$

$$= \sum_{i, j} \alpha_i \alpha_j c_i \delta_{ij} c_j^T$$

$$= \sum_i \alpha_i^2  c_i c_i^T$$

Notice that $D = \sum_i \alpha_i  c_i c_i^T$, $\hat{D} = \sum_i c_i c_i^T$ and $\alpha$ can be 1 or 2 only, we obtain:

$$\qquad DSD = 3D - 2 \hat{D} \ge D$$
