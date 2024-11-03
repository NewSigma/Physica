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
# Exception

## 错误码异常

C++语言使用异常进行错误处理，C语言常使用错误码报告错误。因此在调用C规范的库函数时应使用异常对错误码进行包装。

对任意错误码类型ErrorCodeType，Physica::check使用函数重载进行处理，如有错误则抛出对应的异常对象，其一般实现为:

```
namespace Physica {
    inline void check(ErrorCodeType err) {
        if (err != 0) [[unlikely]]
            throw ExceptionType(err);
    }
}
```
