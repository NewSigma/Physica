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
# Prove of algorithm of HouseHolder

## Theory:
Suppose we have a n-dementional vector $\vec{x} = [x_1, x_2, ..., x_n]^T$ and $\beta = -sgn(x_1) |\vec{x}|$. Define $\vec{a} = \frac{1}{x_1 - \beta} [x_2, x_3, ..., x_n]^T \qquad \tau = 1 - \frac{x_1}{\beta}$.
So we have

$$ (I - \tau \vec{v} \vec{v}^T) \vec{x} = -\beta \hat{e}_1 $$

where $\vec{v}^T = [1, \vec{a}]^T$

Here $sgn(x): \R \rightarrow \R$ is defined as:

$$sgn(x) = \left \{ \begin{matrix} 1 \quad (x \geq 0) \\ -1 \quad (x < 0) \end{matrix} \right.$$

## Prove:

$$ (I - \tau \vec{v} \vec{v}^T) \vec{x} = -\beta \hat{e}_1 $$

$$\Leftrightarrow \vec{x} - \tau \vec{v} \vec{v}^T \vec{x} = -\beta \hat{e}_1$$

Where
$$ \vec{v} \vec{v}^T = \left[ \begin{matrix} 1 & \vec{a}^T \\ \vec{a} & \vec{a} \vec{a}^T \end{matrix} \right]$$

$$ \vec{x} = \left[ \begin{matrix} x_1 \\ (x_1 - \beta) \vec{a} \end{matrix} \right] = \left[ \begin{matrix} x_1 \\ (x_1 + sgn(x_1) |\vec{x}|) \vec{a} \end{matrix} \right] = \left[ \begin{matrix} x_1 \\ sgn(x_1)(|x_1| + |\vec{x}|) \vec{a} \end{matrix} \right] $$

So we have

$$ \vec{v} \vec{v}^T \vec{x} = \left[ \begin{matrix} x_1 + sgn(x_1)(|x_1| + |\vec{x}|) \vec{a}^T \vec{a} \\ (x_1 + sgn(x_1)(|x_1| + |\vec{x}|) \vec{a}^T \vec{a}) \vec{a} \end{matrix} \right] = (x_1 + sgn(x_1)(|x_1| + |\vec{x}|) \vec{a}^T \vec{a}) \vec{v} $$

where

$$ \vec{a}^T \vec{a} = \frac{{x_2}^2 + {x_3}^2 + ... + {x_n}^2}{(x_1 - \beta)^2} = \frac{|\vec{x}|^2 - {x_1}^2}{(|\vec{x}| + |x_1|)^2} = \frac{|\vec{x}| - |x_1|}{|\vec{x}| + |x_1|} $$

$$ \vec{v} \vec{v}^T \vec{x} = (x_1 + sgn(x_1)(|\vec{x}| - |x_1|)) \vec{v} = sgn(x_1) |\vec{x}| \vec{v} $$

$$ \tau = \tau = 1 - \frac{x_1}{\beta} = \frac{|x_1| + |\vec{x}|}{|\vec{x}|} $$

$$ \tau \vec{v} \vec{v}^T \vec{x} = \left[ \begin{matrix} sgn(x_1)(|x_1| + |\vec{x}|) \\ x_2 \\ x_3 \\ ... \\ x_n \end{matrix} \right] $$

So

$$ \vec{x} -  \tau \vec{v} \vec{v}^T \vec{x} = -sgn(x_1)|\vec{x}| \hat{e}_1$$
