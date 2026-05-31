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
# Diff - Automatic differentiation support

## Forward

Implemented using the well-known dual number approach.

## Reverse - Coroutine-based backpropagation

Coroutine-based backpropagation is defined as a backpropagation implementation that uses coroutines to manage the lifetime of the computation graph nodes. *Physica* view forward propagation and backpropagation as a unified process, cross function call boundaries, and expose more optimization opportunities to the compiler.

Physica suspends the coroutine upon completing the forward pass to wait for future gradients, uses RAII to resume coroutine execution upon destruction, and performs gradient accumulation when the coroutine resumes execution.

By [1]:

    Rule  2: What’s good for function values is good for their derivatives.

*Physica* stipulates that if forward propagation does not throw an exception, then backpropagation must also not throw an exception; otherwise, the program is ill-formed.

### Lazy propagation

For types that satisfy concept `ReverseDiff`, the following functions are provided:

``` C++
// Accumulate gradients but don't propagate
void reverse(GradType grad = 1) const noexcept { ... }
```

**Notes**:

1. Although `reverse` is declared `const`, it **does** modify internal state. To see why, consider:

``` C++
VectorND<dfloat> x = ...;
use(x.sum());
```

`sum()` is expected to be `const`. This is natural for forward passes. And the declaration of `reverse` must be consistent. The lowest-level `reverse` implementation has to use `const_cast`.

2. CoDiff:

``` C++
template<class T>
using CoDiff = ...
```

maps type `T` to its corresponding type that satisfies the `ReverseDiff` constraint. Coroutines can only be created on the host side, so the following combinations are illegal:

``` C++
device_obj<CoDiff<...>>
```

### Backpropagation + Expression Templates

In the case of expression templates, passing forward propagation values can avoid the overhead of redundant computations:

``` C++
void reverse(ValueType value, GradType grad = 1) const noexcept { ... }
```

Note that the expression template itself does not perform actual computation nor store the computation results. Therefore, if the result of an expression is needed during the backpropagation process of the expression template, it triggers the necessary forward propagation to obtain that value. This mechanism is equivalent to checkpointing techniques.

## Reference

[1] Evaluating Derivatives (2008); <https://doi.org/10.1137/1.9780898717761>
