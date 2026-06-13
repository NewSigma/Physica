<!--
Copyright 2024-2026 Weibo He.

This file is part of Physica.

Permission is granted to copy, distribute and/or modify this document
under the terms of the GNU Free Documentation License, Version 1.3
or any later version published by the Free Software Foundation;
with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.

You should have received a copy of the GNU Free Documentation License
along with Physica.  If not, see <https://www.gnu.org/licenses/>.
-->
# LinearAlgebra

## Introduction

The *LinearAlgebra* module provides differentiable linear algebra functionalities with SIMD and GPU acceleration.

**Differentiable Scientific Computing**: By extensively and deeply integrating automatic differentiation with modern scientific computing, the goal is to lower the cost of using automatic differentiation techniques and enable continuous benefits from algorithmic advancements through integration with underlying modules. The first step is to have a differentiable linear algebra library.

Inspired by C++ linear algebra libraries such as *Eigen*$^{[1]}$ and *Armadillo*$^{[2]}$, *LinearAlgebra* makes extensive use of expression template techniques to reduce intermediate objects and perform compile-time expression transformations.

## Architecture Design: LValue/RValue Linear Space

Historically, the definitions of lvalue and rvalue in the C language were very simple:

- LValue: An expression that can appear on the left side of an assignment statement. It represents a named memory location with a fixed address.
- RValue: An expression that can only appear on the right side of an assignment statement. It represents a temporary, soon-to-be-destroyed value.

The linear space used for numerical computation is not an abstract linear space; all computations must take place within a specific physical memory. By introducing assignment operations into the linear space, a generalized linear space can be constructed. This design explicitly incorporates expression templates into the type system of *LinearAlgebra*, facilitating the writing of pattern recognition.

The base classes for linear algebra objects generally fall into three categories:

- RValue objects (`RValueVector`, `RValueMatrix`, ...)
- LValue objects (`LValueVector`, `LValueMatrix`, ...)
- Compact objects (`CompactVector`, `CompactMatrix`, ...)

Taking a dense vector as an example, the inheritance hierarchy is: `DenseVector` → `CompactVector` → `LValueVector` → `RValueVector`.

**RValue object**:

Taking an RValue vector as an example, there are only two core operations:

``` C++
Scalar RValueVector::calc(size_t index) { ... }
size_t RValueVector::getLength() { ... }
```

In principle, all "mathematical" vector operations can be implemented on the basis of these core operations. The elements of an rvalue object can be accessed, but their pointers cannot be taken. Therefore, **any** expression template is an rvalue object. The `calc_value()` function is provided to return the value of an element, which is used in cases where gradients are not needed during backpropagation.

**LValue object**:

The only core operation on LValue objects:

``` C++
Scalar* LValueVector::data_ptr(size_t index) { ... }
```

Obviously, LValue objects are calculable, and the operation for calculation is dereferencing:

``` C++
Scalar LValueVector::calc(size_t index) { return *data_ptr(index); }
```

**Compact object**:

The only core operation on compact objects:

``` C++
Scalar* CompactVector::data() { ... }
```

A compact object is one whose elements are continuously distributed in memory:

``` C++
Scalar* CompactVector::data_ptr(size_t index) { return data() + index; }
```

## Concept

Taking a matrix as an example, from the perspective of C language, a matrix is a two-dimensional array:

``` C
void fn(double** matrix) {
    matrix[1][2] = 0;
}
```

Before C++20, from an object-oriented perspective:

``` C++
class Matrix {
    ...
};
```

From C++20 onward, we consider a matrix to be better suited as an abstract *concept*:

``` C++
template<class T>
concept Matrix = ...;
```

We provide four mutually exclusive concepts: `Scalar`, `Vector`, `Matrix`, and `Tensor`. For example: under this design, a `Matrix` with only one column is not considered a `Vector`, which is different from *Eigen*.

## Expression templates

Unary operations can have explicit constraints:

``` C++
template<ExprID, Matrix M>
class UnitaryMatrixExpr { ... };
```

Considering potential non-Abelian properties, binary operations require more general declarations:

``` C++
template<ExprID, class LHS, class RHS>
class BinaryMatrixExpr { ... };
```

Using C++23 explicit object parameter techniques to construct return objects, avoiding a common category of expression template lifetime issues (If you are not familiar with *Eigen*, you are likely to be bitten by this.).

``` C++
auto c = (MatrixND<T>::identity(3) + MatrixND<T>(3, 3, 5)).diag();
// Before C++23: Bad: Heap-use-after-delete
// After C++23: Good: Intermediate results are kept
std::println("{}", c);
```

## Notes

**FastAssign**:

`false`: `assign()` falls back to for-loops on `calc()`/`packet()`
`true`: `assign()` uses optimized algorithm

Motivation: Consider the assignment of linear algebra objects, such as assigning vector:

$$\mathbf{y = x}$$

For dense vectors, element-wise assignment can be done using a for loop, or SIMD can be used with packet as the smallest unit to improve instruction throughput. It is reasonable to assume that a generic implementation using for-loops or SIMD packets provides a predictable performance baseline, bounded by the cost of straightforward element-wise traversal.

There are special cases, including but not limited to sparse structures, symmetry, etc., where template specializations can achieve higher performance than for loops or SIMD implementations. Template specializations are implemented using the `assign` function. `FastAssign = true` indicates that the assign implementation will outperform the for loop or SIMD implementation. `assign` is holistic, and therefore `FastAssign` is propagative.

This option is used to implement heuristic expression transformations.

Example 1 (`assign` outperforms `calc`)
Assigning a sparse vector to a dense vector: a clever implementation uses memset to zero out the dense vector and then fills in the nonzero elements afterward. A trivial implementation does not leverage memset and therefore lacks this optimization.

Example 2 (`calc` is unavailable):
Some calc functions are not implemented for operations that are unlikely to require them. For instance, $\mathbf A^{-1}\mathbf x$ is a vector expression, and one seldom needs just a single element of it. Here, `FastAssign` is true.

**FastPacket**:

`false`: Construction of SIMD objects requires creating a temporary array on the stack  
`true`: No additional overhead  

## Expression template: Why bother?

Q: Handwritten loops can also avoid intermediate object allocations; Why expression templates?  
A: The answer is diverse, one thing for example is: Semantic Preservation Across Call Boundaries  

Consider the implementation of matrix-vector products, a common approach is to lower it to a BLAS implementation (pseudocode!):

``` C++
Vector operator*(Matrix A, Vector x) {
    return BLAS_GEMV(A, x);
}
```

Looks good, right? However, when we actually use it:

``` C++
void somefn(Matrix B, IdentityMatrix I, Vector x) {
    Vector v = (B + I) * x;
    use(v);
}
```

There are a number of special matrices we can use, but here I'm using the identity matrix for clarity. The return value is obviously $B\mathbf x + \mathbf x$, but `operator*` is unaware that matrix $(\mathbf{B + I})$ has special properties, because the function call boundary hinders the propagation of semantics. A BLAS implementation would faithfully iterate over every element of the identity matrix, which is not good. With expression templates, the type `MatrixPlusIdentity` carries the proper semantics:

``` C++
Vector operator*(MatrixPlusIdentity B_plus_I, Vector x) {
    return B_plus_I.B() * x + x;
}
```

where we avoid iterating over the identity matrix and the simplification is zero-overhead. (We really do not want any overhead! Think about $(\mathbf{B + I})^n x$, which appears in exponential matrix functions. If that doesn't scare you, what about $\mathbf{B = C + \lambda I}$ in a shift-and-invert solver?)

Having said that, you do not have to solve the problem using expression templates. But the overall idea is clear: remember to make sufficient use of the information you have!

## Reference

[1] Eigen; <https://eigen.tuxfamily.org>  
[2] Armadillo; <https://arma.sourceforge.net>  
