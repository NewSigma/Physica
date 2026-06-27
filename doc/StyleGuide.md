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
# Physica C++ Style Guide

If not explicitly specified in this document, the Google C++ Style Guide$^{[1]}$ shall be followed by default. For matters not addressed in either this document or [1], consistency with the existing majority of the code shall be maintained.

## C++ standard

No lower than C++23.

## Header file

It is specified that header files have the extension `.h` and source files have the extension `.cpp`.

`<>` are used for standard library and third-party library headers, while `""` are used for *Physica* headers.

### Header guard

All header files should use `#pragma once` to avoid multiple inclusions; the use of `#define` for this purpose is prohibited.

## Scoping

### namespace

In the source file, use `using namespace` to avoid indentation and potential bugs$^{[2]}$

## Class

### Function prototype of `swap`

In practice, it turns out self-swap often leads to bugs. Referring to the `swap` function prototype in [3], it is further specified that `swap` must not perform a self-swap. That is, the general implementation of the `swap` function for any object `T` is as follows:

``` C++
void T::swap(T& obj) noexcept {  
    assert(this != &obj && "[Error]: Self swap is likely a bug");  
    /* Details */  
}
```

## Other C++ features

### RTTI

Prohibited

### `auto`

In function definitions, it is recommended to use trailing return type declarations for return types longer than 4 characters.

## Naming

Template parameters: Encourage the use of abbreviated function templates to eliminate placeholders.

## Formatting

### Namespace Formatting

Indent namespace contents. Avoid nesting namespaces more than 2 levels unless absolutely necessary.

## Reference

[1] Google C++ Style Guide; <https://google.github.io/styleguide/cppguide.html>  
[2] LLVM Coding Standards; <https://llvm.org/docs/CodingStandards.html>  
[3] Scott Meyers 著，侯捷 译. Effective C++：改善程序与设计的55个具体做法[M]. 北京: 电子工业出版社, 2011  
