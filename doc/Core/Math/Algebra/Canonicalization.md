<!--
This file is part of Physica.

SPDX-FileCopyrightText: 2026 Weibo He
SPDX-License-Identifier: GFDL-1.3-no-invariants-or-later
-->

# Canonicalization

Consider, expression $a * b$ for example. Typically, its operands are swappable

$$a * b = b * a,$$

and we call the operation *abelian* mathematically. However, we would like the template parameters of expression templates to follow a particular order for ease of implementation. Here comes the concept of **canonicalization of expression templates**.
