# Diff - 自动微分支持

## Reverse - 协程反向传播

协程反向传播指代使用协程管理计算图生命周期的反向传播实现。C++20的无栈协程为实现具有编译期优化的动态图提供了全新思路。

Physica在前向传播完成时暂停协程以等待未来的梯度。使用RAII在析构时恢复协程执行, 当协程恢复执行时进行梯度累积

由[1]:

    Rule  2: What’s good for function values is good for their derivatives.

规定若前向传播不抛出异常, 则要求反向传播亦不得抛出异常, 否则程序非良构。

提供以下函数

``` C++
void reverse(GradType grad = 1) const noexcept { ... } // 累积梯度但不进行传播
```

对于模板表达式的情况, 传递前向传播的值可以避免重复计算的开销

``` C++
void reverse(ValueType value, GradType grad = 1) const noexcept { ... }
```

## Reference

[1] Evaluating Derivatives (2008); <https://doi.org/10.1137/1.9780898717761>
