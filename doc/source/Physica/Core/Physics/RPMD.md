# Note on normal transformation

In [1], the normal transformation is base on real fourier series, which is unusual in FFT libraries. Here we shall derive complex version of normal transformation.

The hamiltonian reads:

$$H = \sum_{j = 1}^n \frac{p_j^2}{2m} + \frac{1}{2} m \omega^2 (q_{j + 1} - q_j)^2$$

The complex fourier series reads:

$$\tilde{p}_k = \sum_{j = 1}^n p_j \exp(-i \frac{2 \pi}{n}jk)$$

$$p_j = \frac{1}{n} \sum_{k = 0}^{n - 1} \tilde{p}_k \exp(i \frac{2 \pi}{n}jk)$$

Substitude $p_j$ and $q_j$ using $\tilde{p}_k$ and $\tilde{q}_k$, we have

$$n^2 H = \sum_{j \alpha \beta} [\frac{1}{2m} \tilde{p}_\alpha \tilde{p}_\beta + \frac{1}{2} m \omega^2 \tilde{q}_\alpha \tilde{q}_\beta (e^{i \frac{2 \pi}{n} \alpha} - 1) (e^{i \frac{2 \pi}{n} \beta} - 1)] \exp(i \frac{2 \pi}{n} (\alpha + \beta))$$

$$n^2 H = \sum_{\alpha \beta} [\frac{1}{2m} \tilde{p}_\alpha \tilde{p}_\beta + \frac{1}{2} m \omega^2 \tilde{q}_\alpha \tilde{q}_\beta (e^{i \frac{2 \pi}{n} \alpha} - 1) (e^{i \frac{2 \pi}{n} \beta} - 1)] \sum_j \exp(i \frac{2 \pi}{n} (\alpha + \beta))$$

Making use of the fact that $\sum_j \exp(i \frac{2 \pi}{n} (\alpha + \beta)) = \delta_{\alpha + \beta, 0}$, we have

$$n^2 H = \sum_{\alpha} \frac{1}{2m} \tilde{p}_\alpha \tilde{p}_{-\alpha} + \frac{1}{2} m (4 \omega^2 \sin^2(\frac{\pi}{n} \alpha)) \tilde{q}_\alpha \tilde{q}_{-\alpha}$$

$$n^2 H = \sum_{\alpha} \frac{1}{2m} |\tilde{p}_\alpha|^2 + \frac{1}{2} m (4 \omega^2 \sin^2(\frac{\pi}{n} \alpha)) |\tilde{q}_\alpha|^2$$

$$n^2 H = \sum_{\alpha} \frac{1}{2m} |\tilde{p}_\alpha|^2 + \frac{1}{2} m \omega_\alpha^2 |\tilde{q}_\alpha|^2$$

where $\omega_\alpha = 2 \omega \sin(\frac{\pi}{n} \alpha))$

# Reference
[1] M. Ceriotti, M. Parrinello, T. E. Markland and D. E. Manolopoulos, J. Chem. Phys. 133, 124104 (2010).
