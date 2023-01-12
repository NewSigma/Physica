# Note

Learning algorithm $A$ is a random functional that maps training set $X_t$ to a function $f$. The objective of $A$ is minimize the training loss $L_t(X_t; f)$:

$$\mathbb{E} [L_t(X_t,; f)] = \int_{X_t} D[f] P(x) C(x) \mathrm{d}x$$

where $D[f]$ is distance functional that discribe difference between prediction and reality, $P(x)$ is probability that sample x appear, C(x) is cost function that we fail to predict the true value.

Samples $x$ satisfies distribution $G_x$. Given the hyper-param $\lambda$, we have

$$f = A_{\lambda}(X_t) \qquad \forall x \in X_t, x \sim G_x$$

Idea: Training loss $L_t(n)$ is a random process, where $n$ is number of epoch. We assert

$$E[L_t(n_1)] > E[L_t(n_2)]$$
$$\lim_{n \to \infty} D[L_t(n)] = \mathrm{Const}$$

where $n_1 < n_2$.
