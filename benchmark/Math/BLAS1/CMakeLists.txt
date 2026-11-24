# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2025-2026 Weibo He
# SPDX-License-Identifier: MPL-2.0
target_sources(Benchmark PRIVATE Asum.cpp)
target_sources(Benchmark PRIVATE Axpy.cpp)
target_sources(Benchmark PRIVATE Copy.cpp)
target_sources(Benchmark PRIVATE Dot.cpp)
target_sources(Benchmark PRIVATE Nrm2.cpp)
target_sources(Benchmark PRIVATE Scal.cpp)

if(${PHYSICA_CUDA})
    target_sources(Benchmark PRIVATE Axpy.cu)
endif()
