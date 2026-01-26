# Scalar

The Scalar module serves two primary purposes:

- To encapsulate integer and floating-point types within a clean, intuitive interface.
- To enable orthogonal implementation of extended function types - SIMD, complex number and autodiff, with a flexible structure for future additions.

## Notes

- Wrapping floating-point numbers within classes may prevent the compiler from generating Fused Multiply-Add (FMA) instructions. This may or may not align with your expectations. We prefer to avoid aggressive compiler flags such as -ffast-math and instead explicitly use FMA math functions when needed.
- For now, the reverse diff scalar is the only non-copyable scalar.
