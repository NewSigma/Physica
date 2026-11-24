# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2024-2026 Weibo He
# SPDX-License-Identifier: GPL-3.0-or-later
add_subdirectory(BLAS1)
add_subdirectory(BLAS2)
add_subdirectory(BLAS3)
add_subdirectory(Eigen)
add_subdirectory(Matrix)
add_subdirectory(MatrixDecomp)
add_subdirectory(Tensor)
add_subdirectory(VectorExpr)

target_sources(Benchmark PRIVATE
    HouseHolder.cpp
    LnSumExp.cpp
)
