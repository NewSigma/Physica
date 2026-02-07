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
# Exception

在顶层处理复杂数值计算程序的错误通常更具有实用价值。常见数值计算库(LAPACK, cuBLAS)通常使用错误码报告错误，在调用C规范库函数时使用异常进行封装。

对任意错误码类型ErrorCodeType，Physica::check使用函数重载进行处理，如有错误则抛出对应的异常对象，其一般实现为:

``` C++
namespace Physica {
    inline void check(ErrorCodeType err) {
        if (err != 0) [[unlikely]]
            throw ExceptionType(err);
    }
}
```
