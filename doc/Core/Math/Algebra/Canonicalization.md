<!--
Copyright 2026 Weibo He.

This file is part of Physica.

Permission is granted to copy, distribute and/or modify this document
under the terms of the GNU Free Documentation License, Version 1.3
or any later version published by the Free Software Foundation;
with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.

You should have received a copy of the GNU Free Documentation License
along with Physica.  If not, see <https://www.gnu.org/licenses/>.
-->
# Canonicalization

Consider, expression $a * b$ for example. Typically, its operands are swappable

$$a * b = b * a,$$

and we call the operation *abelian* mathematically. However, we would like the template parameters of expression templates to follow a particular order for ease of implementation. Here comes the concept of **canonicalization of expression templates**.
