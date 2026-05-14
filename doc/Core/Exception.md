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

Handling errors at the top level is often more practical for complex numerical computation programs. Common numerical libraries (LAPACK, cuBLAS) typically report errors via error codes, and exceptions are used to wrap C-compatible library function calls.

For any error code type `ErrorCodeType`, `Physica::check` uses function overloading to handle it, throwing the corresponding exception object if an error occurs. The general implementation is:

``` C++
namespace Physica {
    inline void check(ErrorCodeType err) {
        if (err != 0) [[unlikely]]
            throw ExceptionType(err);
    }
}
```
