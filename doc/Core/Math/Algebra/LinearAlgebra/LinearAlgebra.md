<!--
Copyright 2024-2025 Weibo He.

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

受Eigen、Armadillo等C++线性代数库的启发，LinearAlgebra广泛使用模板表达式技术以减少中间对象和进行编译期表达式优化。

## 左值/右值线性空间

历史上, C语言的左值和右值的定义非常简单：

    左值(lvalue)：可以出现在赋值语句左边的表达式。它代表一个有名字、有固定地址的内存位置。
    右值(rvalue)：只能出现在赋值语句右边的表达式。它代表一个临时的、即将消亡的值。

数值计算的线性空间不是抽象的线性空间, 可将赋值运算引入线性空间, 构造一个广义线性空间, 这种构造有助于清晰简洁地用模板表达式实现数学操作。

线性代数对象的基类一般有三种:

右值对象(右值向量RValueVector, 右值矩阵RValueMatrix), 左值对象(LValueVector, LValueMatrix), 连续对象(ContinuousVector, ContinuousMatrix)

以密集向量为例, 继承关系为DenseVector -> ContinuousVector -> LValueVector -> RValueVector

可以计算右值对象的元素但不能取其指针, 任何表达式属于右值对象。

可以取左值对象的元素的指针, 显然左值对象可计算, 计算的操作是解引用

进一步规定连续对象是元素在内存上连续分布的

**右值对象**的核心操作为

1. calc() 计算元素的值, 按值返回
2. getLength(), getRow(), getCol()等 返回线性代数对象的**尺寸**

另外

calc_value() 返回元素的value, 用于反向传播等不需要梯度的情况

## 转换函数

一些常用的逐元素操作被实现为所谓的转换函数，它是返回一个右值对象的成员函数。其尺寸与原对象相同，但在每个元素上进行相同的操作，以右值向量为例:

reals: 返回实部组成的向量  
imags: 虚部  
squaredNorms: 模方  
norms: 模  
values: 值  
grads: 梯度, 对反向传播的纯右值向量调用grads()是不合理的，因为此时计算图尚未生成。调用将导致编译期错误  

## 模板参数

**FastAssign**:

考虑线性代数对象的赋值，如将向量$\mathbf x$赋给向量$\mathbf y$

$$\mathbf{y = x}$$

对于密集向量，使用for逐元素赋值，或使用SIMD以Packet为最小单位以提高指令吞吐量。可以断言，对任意正确的实现，其开销不大于for循环或SIMD实现。

存在特殊情况，包括但不限于稀疏结构、对称性等，使得模板特化的性能高于for循环或SIMD。

模板特化使用函数assign实现，FastAssign = true表示assign实现的性能将优于for或SIMD的实现。

assign具有整体性，因此FastAssign具有传播性。

该选项用于实现启发性表达式变换。

**FastPacket**:

FastPacket = true表示可以快速构造SIMD对象. 对于右值表达式, 需要逐个计算元素并保存在栈上, 此时FastPacket = false

## 模板表达式

一元操作可以显示约束:

``` C++
template<ExprID, Matrix M>
class UnitaryMatrixExpr { ... };
```

考虑到可能的非Abel性, 二元操作需要更一般的声明:

``` C++
template<ExprID, class LHS, class RHS>
class BinaryMatrixExpr { ... };
```

## Reference

[1] Eigen; https://eigen.tuxfamily.org  
[2] Armadillo; https://arma.sourceforge.net  
