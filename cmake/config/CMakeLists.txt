# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2020-2026 Weibo He
# SPDX-License-Identifier: MPL-2.0
include(./PhysicaCheck.cmake)

#CheckRightShift
physica_try_run(CheckRightShift Check/CheckRightShift.cpp)
if(${RUN_RESULT_VAR} MATCHES 1)
    message(STATUS "Performing CheckRightShift --Failed")
else()
    message(STATUS "Performing CheckRightShift --Passed")
endif()

find_package(Hwloc REQUIRED)
physica_cache_probe()

physica_try_compile(RDSEED HostDevAttr/RDSEED.cpp)
string(TOLOWER ${COMPILE_RESULT_VAR} HasRDSEED)

physica_try_compile(RDRAND HostDevAttr/RDRAND.cpp)
string(TOLOWER ${COMPILE_RESULT_VAR} HasRDRAND)

if(${PHYSICA_CUDA})
    physica_cuda_probe(WarpSize)
    physica_cuda_probe(FullUnifiedSupport BOOL)
    physica_cuda_probe(MemoryPoolsSupported BOOL)
endif()

message(STATUS)
file(MAKE_DIRECTORY ${Physica_GENERATED_INCLUDE_DIR}/Physica)
configure_file(Config.h.in ${Physica_GENERATED_INCLUDE_DIR}/Physica/Config.h)
