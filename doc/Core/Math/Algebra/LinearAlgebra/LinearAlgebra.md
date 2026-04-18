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

LinearAlgebra模块提供有SIMD和GPU加速的可微线性代数功能。  

可微分科学计算: 将自动微分广泛而深入的结合现代科学计算, 降低自动微分技术使用成本, 与底层模块结合持续享受算法红利, 首先要有可微分线性代数库。  

受Eigen、Armadillo等C++线性代数库的启发，LinearAlgebra广泛使用表达式模板技术以减少中间对象和进行编译期表达式变换。

## 架构设计: 左值/右值线性空间

历史上, C语言的左值和右值的定义非常简单：

    左值(lvalue)：可以出现在赋值语句左边的表达式。它代表一个有名字、有固定地址的内存位置。
    右值(rvalue)：只能出现在赋值语句右边的表达式。它代表一个临时的、即将消亡的值。

数值计算的线性空间不是抽象的线性空间, 所有计算必须在一定的物理内存中进行。可将赋值运算引入线性空间, 构造广义线性空间。这种设计明确将表达式模板纳入LinearAlgebra的类型系统, 便于模式识别的编写。

线性代数对象的基类一般有三种:

- 右值对象(RValueVector, RValueMatrix, ...)
- 左值对象(LValueVector, LValueMatrix, ...)
- 紧凑对象(CompactVector, CompactMatrix, ...)

以稠密向量为例, 继承关系为DenseVector -> CompactVector -> LValueVector -> RValueVector

**右值对象**:

以右值向量为例, 唯二的核心操作:

``` C++
Scalar RValueVector::calc(size_t index) { ... }
size_t RValueVector::getLength() { ... }
```

原则上所有"数学"向量的操作均可以在核心操作的基础上实现。可以计算右值对象的元素但不能取其指针, 因此**任何**表达式模板都是一个右值对象。提供calc_value()返回元素的value, 用于反向传播中不需要梯度的情况。

**左值对象**:

左值对象的唯一核心操作:

``` C++
Scalar* LValueVector::data_ptr(size_t index) { ... }
```

显然左值对象可计算, 计算的操作是解引用

``` C++
Scalar LValueVector::calc(size_t index) { return *data_ptr(index); }
```

**紧凑对象**:

紧凑对象的唯一核心操作:

``` C++
Scalar* CompactVector::data() { ... }
```

紧凑对象是元素在内存上连续分布的

``` C++
Scalar* CompactVector::data_ptr(size_t index) { return data() + index; }
```

## Concept

以矩阵为例, 在C语言观点下, 矩阵是一个二维数组:

``` C
void fn(double** matrix) {
    matrix[1][2] = 0;
}
```

C++20以前, 面向对象的观点:

``` C++
class Matrix {
    ...
};
```

C++20及以后, 我们认为矩阵更适合作为抽象的concept(概念):

``` C++
template<class T>
concept Matrix = ...;
```

我们提供四个互斥的concept: `Scalar`, `Vector`, `Matrix`, `Tensor`; 例: 这个设计下, 只有一列的`Matrix`不是一个`Vector`。

## 表达式模板

一元操作可以显式约束:

``` C++
template<ExprID, Matrix M>
class UnitaryMatrixExpr { ... };
```

考虑到可能的非Abel性, 二元操作需要更一般的声明:

``` C++
template<ExprID, class LHS, class RHS>
class BinaryMatrixExpr { ... };
```

使用C++23显式对象参数技术构造返回对象, 避免一类表达式模板常见的生命周期问题:

``` C++
auto c = (MatrixND<T>::identity(3) + MatrixND<T>(3, 3, 5)).diag();
// Before C++23: Bad: Heap-use-after-delete
// After C++23: Good: Intermediate results are kept
std::println("{}", c);
```

## Notes

**FastAssign**:

考虑线性代数对象的赋值，如将向量$\mathbf x$赋给向量$\mathbf y$

$$\mathbf{y = x}$$

对于稠密向量，使用for逐元素赋值，或使用SIMD以Packet为最小单位以提高指令吞吐量。可以断言，对任意正确的实现，其开销不大于for循环或SIMD实现。

存在特殊情况，包括但不限于稀疏结构、对称性等，使得模板特化的性能高于for循环或SIMD。

模板特化使用函数assign实现，FastAssign = true表示assign实现的性能将优于for或SIMD的实现。

assign具有整体性，因此FastAssign具有传播性。

该选项用于实现启发性表达式变换。

**FastPacket**:

`false`: SIMD对象的构造需要在栈上创建临时数组
`true`: 无额外开销

## Reference

[1] Eigen; <https://eigen.tuxfamily.org>  
[2] Armadillo; <https://arma.sourceforge.net>  
