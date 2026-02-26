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
# Diff - 自动微分支持

## Forward

使用常规的对偶数方法实现

## Reverse - 协程反向传播

协程反向传播指代使用协程管理计算图生命周期的反向传播实现。我们将前向传播和反向传播视为一个完整的过程, 跨越函数调用边界, 为编译器暴露更多优化机会。

Physica在前向传播完成时暂停协程以等待未来梯度, 使用RAII在析构时恢复协程执行, 当协程恢复执行时进行梯度累积。

由[1]:

    Rule  2: What’s good for function values is good for their derivatives.

规定若前向传播不抛出异常, 则要求反向传播亦不得抛出异常, 否则程序非良构。

对满足ReferseDiff概念的类型, 提供以下函数

    ``` C++
    // 累积梯度但不进行传播
    void reverse(GradType grad = 1) const noexcept { ... }
    ```

考虑

    ``` C++
    VectorND<dfloat> x = ...;
    use(x.sum());
    ```

我们希望将`sum()`声明为`const`, `reverse`的声明应当保持一致:

    ``` C++
    // const reverse: working
    // non-const reverse: does not compile
    x.sum().reverse();
    ```

最底层的`reverse`实现需要使用`const_cast`.

### CoDiff

    ``` C++
    tempalate<class T>
    using CoDiff = ...
    ```

将类型 T 映射为其满足 ReverseDiff 约束的对应类型。协程只能在主机端(host)创建，因此下列组合属于非法情况：

    ``` C++
    device_obj<CoDiff<...>>
    ```

### 反向传播 + 表达式模板

对于表达式模板的情况, 传递前向传播的值可以避免重复计算的开销

    ``` C++
    void reverse(ValueType value, GradType grad = 1) const noexcept { ... }
    ```

注意到，表达式模板本身并不执行实际计算, 亦不储存计算结果。因此，表达式模板反向传播过程中若需使用表达式的结果，会触发必要的前向传播以获取该值，该机制等同于Checkpoint技术。

## Reference

[1] Evaluating Derivatives (2008); <https://doi.org/10.1137/1.9780898717761>
