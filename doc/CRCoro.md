# 奇异递归协程

## Motivation

C++20无栈协程目前已经相对成熟，但是要将旧代码库迁移到新的语言标准通常会遇到一些困难。具体到模板的迁移，我们可能尝试写出这样的代码

```C++
template<bool Coro>
Result fn() {
    if constexpr (Coro)
        co_return foo();
    else
        return bar();
}
```

遗憾的是，上述代码不能通过编译。因为C++标准规定`return`关键字不能和协程相关关键词`co_return`, `co_yield`, `co_await`出现在同一函数中。对于这种函数的迁移，我们只能采取妥协的办法，本文提出一种结合奇异递归模板模式和协程的非侵入的范式尝试解决该问题。

## Implementation

奇异递归协程的核心设计思想是要求普通函数无开销或小开销地转化为一个平庸的协程(依赖于编译器的优化)，最小代码如下

```C++
#include <utility>
#include <coroutine>

template<class T>
struct CRCoro { // Curiously Recurring Coroutine
    struct RValueWrapper {
        T* p;

        operator T&&() const noexcept { return std::move(*p); }
    };

    using promise_type = T;

    T& getDerived() noexcept { return *static_cast<T*>(this); }

    auto get_return_object() noexcept { return RValueWrapper(&getDerived()); }
    std::suspend_never initial_suspend() noexcept { return {}; }
    std::suspend_never final_suspend() noexcept { return {}; }
    void return_value(T&& x) noexcept {
        getDerived() = std::move(x);
    }
    void unhandled_exception() {}
};
```

要设计一个协程，我们需要规定以下函数:

`get_return_object`: 返回对象，必选
`initial_suspend`: 初始暂停行为，必选
`final_suspend`: 最终暂停行为，必选
`return_value`/`return_void`: 返回值/返回空，任选其一
`yield_value`: 产出值，可选
`unhandled_exception`: 发生异常时的行为，必选

一个平庸的协程不应包含任何等待点, 因此我们要求`initial_suspend`和`final_suspend`返回`std::suspend_never`, 不实现`yield_value`。
我们希望协程返回协程体产生的值, 因此在return的两个版本中选择`return_value`, 返回值必须保存在承诺对象中。
在返回对象的实现中, 我们构造了右值包装对象`RValueWrapper`, 通过移动构造的方式将承诺对象返回给调用者。
我们选择不处理异常, `unhandled_exception`为空

具体使用时，我们要求返回对象继承自`CRCoro`，正如我们习惯在奇异递归模板模式中做的事情:

```C++
#include <cstdio>

struct A : public CRCoro<A> {
    int a;
};

A fn() {
    A aa{};
    aa.a = 5;
    co_return std::move(aa);
}

int main() {
    printf("%d\n", fn().a); // Expect printing 5
    return 0;
}
```

不幸的是，上述代码含有一个UB。三大主要编译器中只有MSVC能够给出预期的结果，GCC和Clang都会选择在协程栈帧析构后进行返回对象的转换$^{[1]}$。这一UB有希望在CWG2563$^{[2]}$中修复。在正式修复以前，我们可以采用以下实现规避该UB:

```C++
#include <utility>
#include <coroutine>

template<class T>
struct CRCoro { // Curiously Recurring Coroutine
    struct RValueWrapper {
        T* p;

        ~RValueWrapper() {
            // Workaround, wait for CWG2563
            std::coroutine_handle<T>::from_promise(*p).destroy();
        }

        operator T&&() const noexcept { return std::move(*p); }
    };

    using promise_type = T;

    T& getDerived() noexcept { return *static_cast<T*>(this); }

    auto get_return_object() noexcept { return RValueWrapper(&getDerived()); }
    std::suspend_never initial_suspend() noexcept { return {}; }
    std::suspend_always final_suspend() noexcept { return {}; }
    void return_value(T&& x) noexcept {
        getDerived() = std::move(x);
    }
    void unhandled_exception() {}
};
```

注意该Workaround将会阻碍Clang的协程消除优化。

## Conclusion

我们提出奇异递归协程帮助迁移模板代码。一个从不等待的协程是非常"奇异"的, 它也未必是解决迁移问题的最佳方案。从另一角度, 我们提供了一个有趣的corner case。这对这一case，我们可以问两个问题

1. C++标准如何定义该case的行为
2. 编译器如何优化该case

## Reference

[1] 118074; https://gcc.gnu.org/bugzilla/show_bug.cgi?id=118074
[2] CWG2563; https://cplusplus.github.io/CWG/issues/2563.html
