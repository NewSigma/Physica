# Some notes on function splitOffTwoRows()

Given a 2x2 matrix:
$$ A = \left[ \begin{matrix} a & b \\ c & d \end{matrix} \right] ,$$
whose eigenvalues are $ \frac{a + d}{2}  \pm \sqrt{\frac{(a - d)^2}{4} + bc}$

Now we want a orthogonal matrix $U$, such that

$$UAU^T = M$$

where $M$ is a upper triangular matrix.

Suppose $U$ is a givens matrix, that is:

$$ U = \left[ \begin{matrix} \cos{\theta} & \sin{\theta} \\ -\sin{\theta} & \cos{\theta} \end{matrix} \right] $$

So the (2, 1) element of $M$ is

$$ M(2, 1) = c \cos^2{\theta} - b \sin^2{\theta} + (d - a) \sin{\theta} \cos{\theta} $$

Let $M(2, 1) = 0$, note that $\sin{\theta}$ cannot be 0, or matrix $U$ will be identity, we have

$$ c \cot^2{\theta} + (d - a) \cot{\theta} - b = 0 $$

Solve the equation, we have
$$\cot{\theta} = \frac{\frac{a - d}{2}  \pm \sqrt{\frac{(a - d)^2}{4} + bc}}{c}$$

which is similiar to the expression of eigenvalues.
