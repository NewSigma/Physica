<!--
Copyright 2024 Weibo He.

This file is part of Physica.

Permission is granted to copy, distribute and/or modify this document
under the terms of the GNU Free Documentation License, Version 1.3
or any later version published by the Free Software Foundation;
with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.

You should have received a copy of the GNU Free Documentation License
along with Physica.  If not, see <https://www.gnu.org/licenses/>.
-->
# Physica C++编程规范

如本规范未明确指出，默认遵循Google C++编程规范$^{[1]}$。本规范与[1]均未指出的，应与已有的多数代码一致。

## C++版本

目前使用的C++版本应不低于C++20。

## 头文件

规定头文件拓展名为.h，源文件拓展名为.cpp

标准库、第三方库头文件使用<>, Physica头文件使用""

### Header guard

所有头文件应使用#pragma once以避免多重include，禁止使用#define

## Scoping

### namespace

源文件中使用using namespace以避免缩进

## Class

### swap函数原型

参考[2]中的swap函数原型，进一步规定swap不可与自身发生交换，即任意对象T的swap函数的一般实现为

``` C++
void T::swap(T& __restrict obj) noexcept {  
    assert(this != &obj && "[Error]: Self swap is likely a bug");  
    /* Details */  
}
```

实践中发现Self swap常导致Bug，通过添加__restrict帮助编译器优化和提供编译期Self swap警告

## 其他C++特性

### 静态断言

类模板的静态断言应尽可能放置在Traits中

### RTTI

不得使用

## 命名

Namespace: 大驼峰式

## 单元测试

单元测试用时一般不超过1 min，合理情况下可适当延长，任何情况下都不应超过10 min。

## Reference

[1] Google C++ Style Guide; https://google.github.io/styleguide/cppguide.html  
[2] Scott Meyers 著，侯捷 译. Effective C++：改善程序与设计的55个具体做法[M]. 北京: 电子工业出版社, 2011  
