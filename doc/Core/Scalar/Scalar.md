<!--
Copyright 2025-2026 Weibo He.

This file is part of Physica.

Permission is granted to copy, distribute and/or modify this document
under the terms of the GNU Free Documentation License, Version 1.3
or any later version published by the Free Software Foundation;
with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.

You should have received a copy of the GNU Free Documentation License
along with Physica.  If not, see <https://www.gnu.org/licenses/>.
-->
# Scalar

The Scalar module serves two primary purposes:

- To encapsulate integer and floating-point types within a uniform interface.
- To enable orthogonal implementation of extended function types - SIMD, complex number and autodiff, with a flexible structure for future additions.

## Notes

- Wrapping floating-point numbers within classes may prevent the compiler from generating Fused Multiply-Add (FMA) instructions. This may or may not align with your expectations. We prefer to avoid aggressive compiler flags such as -ffast-math and instead explicitly use FMA math functions when needed.
- For now, the reverse diff scalar is the only non-copyable scalar.
