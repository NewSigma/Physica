# Notes on HFSolver

In a general HF procedure, we require the wave function is normalized:

$$c^T S c = 1$$

where $c$ is the wave function vector, $S$ is the overlap matrix.

The passage aim to clarify the normalization is unnecessary in the code.

## Prove

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
