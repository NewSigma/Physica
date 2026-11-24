# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2026 Weibo He
# SPDX-License-Identifier: MPL-2.0
add_subdirectory(Executor)
add_subdirectory(Future)
add_subdirectory(MPI)

target_sources(PhysicaCore PRIVATE
    Task/Task.cpp
    MPI.cpp
    SubProcess.cpp
    ThreadPool.cpp
)

if(${PHYSICA_CUDA})
    target_sources(PhysicaCore PRIVATE
        CUDAContext.cu
        CUDAEvent.cu
        CUDAStream.cu
        KernelConfig.cu
    )
endif()
