# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2026 Weibo He
# SPDX-License-Identifier: GPL-3.0-or-later
target_sources(Benchmark PRIVATE GEMV.cpp)

if(${PHYSICA_CUDA})
    target_sources(Benchmark PRIVATE GEMV.cu)
endif()
