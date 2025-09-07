<!--
Copyright 2024-2025 Weibo He.

This file is part of Physica.

Permission is granted to copy, distribute and/or modify this document
under the terms of the GNU Free Documentation License, Version 1.3
or any later version published by the Free Software Foundation;
with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.

You should have received a copy of the GNU Free Documentation License
along with Physica.  If not, see <https://www.gnu.org/licenses/>.
-->
# Algorithm of HouseHolder

We are following a slightly different conversion from that of textbook algorithm$^{[1]}$.

## Theory

Let n-dementional vector $\mathbf{x} = [x_1, x_2, ..., x_n]^T$. We have

$$(I - \tau \mathbf{v} \mathbf{v}^T) \mathbf{x} = -\text{sgn}(x_1) |\mathbf{\mathbf{x}}| \hat{e}_1,$$

where $\mathbf{v}^T = [1, \mathbf{a}]^T$, $\mathbf{\mathbf{a}} = \frac{1}{x_1 + \text{sgn}(x_1) |\mathbf{\mathbf{x}}|} [x_2, x_3, ..., x_n]^T$, $\tau = 1 + |x_1|/|\mathbf{x}|$. Sign function is defined as

$$\text{sgn}(x) = \left\{ \begin{matrix} 1 \quad (x \geq 0) \\ -1 \quad (x < 0) \end{matrix} \right.$$

Complex version is obtained simply by replacing $\text{sgn}(x)$ with phase $x/|x|$.

## Prove

$$(I - \tau \mathbf{v} \mathbf{v}^T) \mathbf{x} = -\text{sgn}(x_1) |\mathbf{\mathbf{x}}| \hat{e}_1$$

$$\Rightarrow \mathbf{x} - \tau \mathbf{v} \mathbf{v}^T \mathbf{x} = -\text{sgn}(x_1) |\mathbf{\mathbf{x}}| \hat{e}_1,$$

where

$$\mathbf{v} \mathbf{v}^T = \left[ \begin{matrix} 1 & \mathbf{a}^T \\ \mathbf{a} & \mathbf{a} \mathbf{a}^T \end{matrix} \right],$$

$$\mathbf{x} = \left[ \begin{matrix} x_1 \\ (x_1 + \text{sgn}(x_1) |\mathbf{\mathbf{x}}|) \mathbf{a} \end{matrix} \right] = \left[ \begin{matrix} x_1 \\ (x_1 + \text{sgn}(x_1) |\mathbf{x}|) \mathbf{a} \end{matrix} \right] = \left[ \begin{matrix} x_1 \\ \text{sgn}(x_1)(|x_1| + |\mathbf{x}|) \mathbf{a} \end{matrix} \right].$$

Substitude into the above formula

$$\mathbf{v} \mathbf{v}^T \mathbf{x} = \left[ \begin{matrix} x_1 + \text{sgn}(x_1)(|x_1| + |\mathbf{x}|) \mathbf{a}^T \mathbf{a} \\ (x_1 + \text{sgn}(x_1)(|x_1| + |\mathbf{x}|) \mathbf{a}^T \mathbf{a}) \mathbf{a} \end{matrix} \right] = (x_1 + \text{sgn}(x_1)(|x_1| + |\mathbf{x}|) \mathbf{a}^T \mathbf{a}) \mathbf{v},$$

where

$$\mathbf{a}^T \mathbf{a} = \frac{{x_2}^2 + {x_3}^2 + ... + {x_n}^2}{(x_1 \text{sgn}(x_1) |\mathbf{\mathbf{x}}|)^2} = \frac{|\mathbf{x}|^2 - {x_1}^2}{(|\mathbf{x}| + |x_1|)^2} = \frac{|\mathbf{x}| - |x_1|}{|\mathbf{x}| + |x_1|}.$$

As a result

$$\mathbf{v} \mathbf{v}^T \mathbf{x} = (x_1 + \text{sgn}(x_1)(|\mathbf{x}| - |x_1|)) \mathbf{v} = \text{sgn}(x_1) |\mathbf{x}| \mathbf{v},$$

and

$$\tau \mathbf{v} \mathbf{v}^T \mathbf{x} = \left[ \begin{matrix} x_1 + \text{sgn}(x_1)|\mathbf{x}| \\ x_2 \\ x_3 \\ ... \\ x_n \end{matrix} \right].$$

Finally

$$\mathbf{x} -  \tau \mathbf{v} \mathbf{v}^T \mathbf{x} = -\text{sgn}(x_1)|\mathbf{x}| \hat{e}_1$$

## References

[1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013:236
