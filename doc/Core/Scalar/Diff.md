# Diff

## Reverse - 协程反向传播

Reference [1]:

    Rule  2: What’s good for function values is good for their derivatives.

若前向传播不抛出异常, 则要求反向传播亦不抛出异常, 否则程序非良构. 使用RAII在析构时进行反向传播, 使用协程管理计算图生命周期, 并提供以下函数

    void reverse(GradType grad = 1) const noexcept { ... } // 累计梯度但不进行传播

对于模板表达式的情况, 传递前向传播的值可以避免重复计算的开销

    void reverse(ValueType value, GradType grad = 1) const noexcept { ... }

## Reference

[1] Evaluating Derivatives (2008); https://doi.org/10.1137/1.9780898717761
