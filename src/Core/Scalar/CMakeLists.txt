# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2025-2026 Weibo He
# SPDX-License-Identifier: GPL-3.0-or-later
add_subdirectory(IntegerImpl)
add_subdirectory(RealImpl)

target_sources(PhysicaCore PRIVATE
    ArrayArithmetic.cpp
    Integer.cpp
    Rational.cpp
    Real.cpp
)
