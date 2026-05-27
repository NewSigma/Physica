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
# Curiously Recurring Coroutine

## Motivation

C++20 stackless coroutines are currently relatively mature, but migrating legacy codebases to the new language standard often comes with some difficulties. Specifically regarding the migration of templates, we might try to write code like this:

```C++
template<bool Coro>
Result fn() {
    if constexpr (Coro)
        co_return foo();
    else
        return bar();
}
```

Unfortunately, the above code does not compile. This is because the C++ standard prohibits the `return` keyword and coroutine keywords from appearing in the same function. For migrating such functions, we have no choice but to adopt a compromise — there are simply too many of them. This paper proposes a non-intrusive paradigm that combines the Curiously Recurring Template Pattern (CRTP) and coroutines in an attempt to solve this problem.

## Implementation

The core design philosophy of *curiously recurring coroutine* is to require ordinary functions to be transformed into a trivial coroutine with zero or low overhead (relying on compiler optimizations), with minimal code as follows:

```C++
#include <utility>
#include <coroutine>

template<class T>
struct CRCoro { // Curiously Recurring Coroutine
    using promise_type = T;

    T& getDerived() noexcept { return *static_cast<T*>(this); }

    auto get_return_object() noexcept {
        struct RValueWrapper {
            This* p;

            operator T&&() const noexcept { return std::move(p->getDerived()); }
        };
        return RValueWrapper(this);
    }
    void await_transform(auto&&) noexcept = delete;
    std::suspend_never initial_suspend() noexcept { return {}; }
    std::suspend_never final_suspend() noexcept { return {}; }
    void return_value(T&& x) noexcept { getDerived() = std::move(x); }
    void unhandled_exception() {}
};
```

To design a coroutine, we need to specify the following functions:

`get_return_object`: returns an object, required  
`initial_suspend`: initial suspension behavior, required  
`final_suspend`: final suspension behavior, required  
`return_value/return_void`: return value / return void, choose one  
`yield_value`: yield value, optional  
`unhandled_exception`: unhandled exception, required  

A trivial coroutine should not contain any suspension points, so we require that `initial_suspend` and `final_suspend` return `std::suspend_never`, and that `yield_value` is not implemented.  
We want the coroutine to return the value produced by the coroutine body, so we choose `return_value` from the two return versions, and the return value must be stored in the promise object.  
In the implementation of `get_return_object`, we construct a wrapper object `RValueWrapper`, which returns the promise object to the caller via move construction.  
We choose not to handle exceptions, so `unhandled_exception` is left empty.  

In specific usage, we require the coroutine return type `A` to inherit from `CRCoro<A>`, just as we are accustomed to doing in CRTP:

```C++
#include <cassert>
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
    assert(fn().a == 5);
    return 0;
}
```
