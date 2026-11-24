# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2026 Weibo He
# SPDX-License-Identifier: GPL-3.0-or-later
add_subdirectory(CUDA)
add_subdirectory(Unix)

if (${PHYSICA_MIMALLOC})
    target_sources(PhysicaCore PRIVATE Allocator/Mimalloc.cpp)
endif()

target_sources(PhysicaCore PRIVATE
    BitArray.cpp
    Cycler.cpp
    Memory.cpp
    NoImpl.cpp
    RingBuffer.cpp
    Suspend.cpp
)
