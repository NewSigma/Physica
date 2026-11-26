<!--
This file is part of Physica.

SPDX-FileCopyrightText: 2020-2026 Weibo He
SPDX-License-Identifier: MPL-2.0
-->

# Physica

[![License](https://img.shields.io/badge/License-MPL_2.0-blue.svg)](https://opensource.org/license/mpl-2-0) [![Linux](https://img.shields.io/badge/os-Linux-green.svg)](https://shields.io/)

*Physica* is a high-performance and scalable C++ template library, dedicated to providing a rapid-iteration platform for new algorithms.

It is inspired by [*Eigen*](https://eigen.tuxfamily.org) but integrates differentiation, tensor expressions and GPU execution into a uniform type system. This makes it elegant to write a numerical algorithm, differentiate it, and parallelize it.

## Why Physica

- **Multilinear algebra**: `Scalar`, `Vector`, `Matrix` and `Tensor` are distinct concepts with expression support.
- **Expression templates**: `A * x + b` remains a lazy expression, enabling operator fusion and compile-time expression transformations.
- **Automatic differentiation**: forward- and reverse-mode automatic differentiation are available for scalars, vectors, matrices and tensors, not as an external layer.
- **GPU execution**: host and device objects share the same interface via `device_obj<T>` and `toDevice()`.
- **Lifetime safety**: expression nodes own their rvalue operands, avoiding the dangling pitfalls of traditional expression templates.
- **Native scientific computing**: ODE/PDE, statistics, optimization and computational physics components build on the same foundation.

## Quickstart

The following example computes `ReLU(Ax + b)` and propagates the gradient of `loss` back to `A`, `x` and `b`:

```cpp
using T = float64;
using dfloat = Diff<T, DiffMode::Reverse>;

int main() {
    MatrixND<dfloat> A{{1, 2}, {-1, 1}};
    VectorND<dfloat> x{1, 1};
    VectorND<dfloat> b{-1, -0.5};

    auto y = relu(A * x + b);
    auto loss = square(y).sum();
    loss.reverse_final();

    std::println("dL/dA =\n{}", A.grads());
    std::println("dL/dx = {}", x.grads());
    std::println("dL/db = {}", b.grads());
}
```

The complete example is available at [`examples/Core/Diff.cpp`](examples/Core/Diff.cpp).

## Usage

Refer to [Install.md](doc/Install.md) for a detailed guide on using *Physica*.

## Maintainers

Weibo He (<NewSigma@163.com>)  

## License

[MPL-2.0](LICENSE) © Weibo He
