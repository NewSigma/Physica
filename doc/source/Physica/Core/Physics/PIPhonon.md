# How to apply translation invariance

The translation invariance requires that

$$\sum_i F_{i \sigma} = 0 \qquad \sum_i p_{i \sigma} = 0$$,

where $i$ is index of atoms in unit cell and $\sigma = x, y, z$.

Take the momentum correlation matrix $M$ as example, which is defined as:

$$\langle i \sigma|M|j \sigma' \rangle = \sum_k p^k_{i \sigma} p^k_{j \sigma'}$$

Apply the translation invariance:

$$\sum_i (p_{i \sigma}^k + \Delta p_\sigma^k) = 0$$

Our goal is deltas of matrix elements

$$\langle i \sigma|\Delta M|j \sigma' \rangle = \sum_k(p^k_{i \sigma} \Delta p^k_{\sigma'} + p^k_{j \sigma'} \Delta p^k_{\sigma} + \Delta p^k_{\sigma} \Delta p^k_{\sigma'})$$

$$= \sum_k[p^k_{i \sigma} (-\frac{1}{n} \sum_\alpha p^k_{\alpha \sigma'}) + p^k_{j \sigma'} (-\frac{1}{n} \sum_\alpha p^k_{\alpha \sigma}) + (-\frac{1}{n} \sum_\alpha p^k_{\alpha \sigma}) (-\frac{1}{n} \sum_\alpha p^k_{\alpha \sigma'})]$$

$$= \frac{1}{n^2} \sum_{\alpha \beta} \sum_k p^k_{\alpha \sigma} p^k_{\beta \sigma'} - \frac{1}{n} \sum_\alpha \sum_k p^k_{i \sigma} p^k_{\alpha \sigma'} - \frac{1}{n} \sum_\alpha \sum_k p^k_{\alpha \sigma} p^k_{j \sigma'}$$

$$= \frac{1}{n^2} \sum_{\alpha \beta} \langle \alpha \sigma|M|\beta \sigma' \rangle - \frac{1}{n} \sum_{\alpha} \langle i \sigma|M|\alpha \sigma' \rangle - \frac{1}{n} \sum_{\alpha} \langle \alpha \sigma|M|j \sigma' \rangle$$
