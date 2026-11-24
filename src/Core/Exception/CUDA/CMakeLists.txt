# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2026 Weibo He
# SPDX-License-Identifier: MPL-2.0
target_sources(PhysicaCore PRIVATE
    cuBLAS.cu
    CUDA.cu
    cuDSS.cu
    cuRAND.cu
    cuSolver.cu
)
