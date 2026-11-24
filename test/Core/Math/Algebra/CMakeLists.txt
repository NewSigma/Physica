# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2024-2026 Weibo He
# SPDX-License-Identifier: GPL-3.0-or-later
add_subdirectory(LinearAlgebra)
physica_add_test(Canonicalization Canonicalization.cpp)
physica_add_test(Polynomial Polynomial.cpp)
