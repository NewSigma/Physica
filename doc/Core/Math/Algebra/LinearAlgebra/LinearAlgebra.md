# LinearAlgebra

## Introduction

受Eigen、Armadillo等C++线性代数库的启发，Physica线性代数模块广泛使用模板表达式技术以减少中间对象和进行编译期表达式优化。

## 模板参数

**FastAssign**

考虑线性代数对象的赋值，如将向量$\mathbf x$赋给向量$\mathbf y$

$$\mathbf{y = x}$$

对于密集向量，使用for逐元素赋值，或使用SIMD以Packet为最小单位以提高指令吞吐量。可以断言，对任意正确的实现，其开销不大于for循环或SIMD实现。

存在特殊情况，包括但不限于稀疏结构、对称性等，使得模板特化的性能高于for循环或SIMD。

模板特化使用函数assignTo实现，FastAssign = true表示assignTo实现的性能将优于for或SIMD的实现。

assignTo具有整体性，因此FastAssign具有传播性。

该选项可用于实现启发性表达式变换。

## Reference

[1] Eigen; https://eigen.tuxfamily.org/
[2] Armadillo; https://arma.sourceforge.net/
