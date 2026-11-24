# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2026 Weibo He
# SPDX-License-Identifier: MPL-2.0
target_sources(PhysicaCore PRIVATE ProcessFuture.cpp)

if(${PHYSICA_CUDA})
    target_sources(PhysicaCore PRIVATE StreamFuture.cu)
endif()
