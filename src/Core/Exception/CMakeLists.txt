# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2025-2026 Weibo He
# SPDX-License-Identifier: GPL-3.0-or-later
target_sources(PhysicaCore PRIVATE Syscall.cpp SystemError.cpp)

if (${PHYSICA_CUDA})
    add_subdirectory(CUDA)
endif()

if (${PHYSICA_MKL})
    add_subdirectory(MKL)
endif()

if (${PHYSICA_MPI})
    target_sources(PhysicaCore PRIVATE MPI.cpp)
endif()
