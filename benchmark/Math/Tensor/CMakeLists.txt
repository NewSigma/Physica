# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2026 Weibo He
# SPDX-License-Identifier: MPL-2.0
if(${PHYSICA_CUDA})
    target_sources(Benchmark PRIVATE StridedTensor.cu)
endif()
