<!--
This file is part of Physica.

SPDX-FileCopyrightText: 2024 Weibo He
SPDX-License-Identifier: GFDL-1.3-no-invariants-or-later
-->

# ProbDistribution

ProbDistribution computes the probability distribution function of a one-dimensional random variable using sampling semantics. Consider the random variable:

$$y = \ln\left( 1 + \sqrt{1-e^{-1/n}} \right) \sum_{i = 1}^n s_i$$

where the independent random variables $s_i = \pm 1$ follow a binomial distribution with $p = 1/2$.

![](./ProbDistribution.png)

**Fig. 1** Blue: probability density distribution of the random variable $y$; Green: standard normal distribution

All moments of the random variable $y$ equal those of a standard normal variable, and its moment-generating function exists, so $y \sim N(0, 1)$ as $n \to \infty$.
