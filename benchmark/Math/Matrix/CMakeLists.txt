# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2026 Weibo He
# SPDX-License-Identifier: MPL-2.0
target_sources(Benchmark PRIVATE IdentityMV.cpp)
target_sources(Benchmark PRIVATE MatrixPow.cpp)
target_sources(Benchmark PRIVATE MatrixTrig.cpp)
target_sources(Benchmark PRIVATE CompactMatrix.cpp)
target_sources(Benchmark PRIVATE SyMV.cpp)
target_sources(Benchmark PRIVATE Transpose.cpp)

if(${PHYSICA_CUDA})
    target_sources(Benchmark PRIVATE StridedMatrix.cu)
endif()
