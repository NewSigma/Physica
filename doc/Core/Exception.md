<!--
This file is part of Physica.

SPDX-FileCopyrightText: 2024-2026 Weibo He
SPDX-License-Identifier: GFDL-1.3-no-invariants-or-later
-->

# Exception

Handling errors at the top level is often more practical for complex numerical computation programs. Common numerical libraries (LAPACK, cuBLAS) typically report errors via error codes, and exceptions are used to wrap C-compatible library function calls.

For any error code type `ErrorCodeType`, `Physica::check` uses function overloading to handle it, throwing the corresponding exception object if an error occurs. The general implementation is:

``` C++
void Physica::check_backend(int err) {
    constinit static const ErrorCategoryImpl instance{};
    if (err != 0) [[unlikely]]
        throw SystemError(err, instance);
}
```
