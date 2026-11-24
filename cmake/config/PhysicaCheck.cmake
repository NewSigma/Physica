# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2020-2026 Weibo He
# SPDX-License-Identifier: MPL-2.0
if(${PHYSICA_CUDA})
    set(PHYSICA_CHECK_LIBS CUDA::cudart)
endif()

macro(physica_try_compile NAME SOURCE)
    try_compile(COMPILE_RESULT_VAR
                SOURCES ${CMAKE_HOME_DIRECTORY}/cmake/config/${SOURCE}
                CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${CUDAToolkit_INCLUDE_DIRS};${HWLOC_INCLUDE_DIRS}" -DCOMPILE_DEFINITIONS=-march=${PHYSICA_ARCH}
                LINK_LIBRARIES ${PHYSICA_CHECK_LIBS} ${HWLOC_LIBRARIES}
                OUTPUT_VARIABLE COMPILE_LOG)
endmacro()

macro(physica_try_run NAME SOURCE)
    try_run(RUN_RESULT_VAR COMPILE_RESULT_VAR ${CMAKE_BINARY_DIR} ${CMAKE_HOME_DIRECTORY}/cmake/config/${SOURCE}
            CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${CUDAToolkit_INCLUDE_DIRS};${HWLOC_INCLUDE_DIRS}" -DCOMPILE_DEFINITIONS=-march=${PHYSICA_ARCH}
            LINK_LIBRARIES ${PHYSICA_CHECK_LIBS} ${HWLOC_LIBRARIES}
            COMPILE_OUTPUT_VARIABLE COMPILE_LOG
            RUN_OUTPUT_STDOUT_VARIABLE RUN_OUTPUT_VAR)

    if(NOT ${COMPILE_RESULT_VAR})
        message(FATAL_ERROR "Failed to compile ${NAME}: ${COMPILE_LOG}")
    endif()

    if(${RUN_RESULT_VAR} MATCHES FAILED_TO_RUN)
            message(FATAL_ERROR "Failed to run ${NAME}")
    endif()
endmacro()

macro(physica_cuda_probe NAME)
    if(NOT DEFINED ${NAME})
        physica_try_run(${NAME} CUDADevAttr/${NAME}.cpp)
        if("${ARGV1}" STREQUAL "BOOL")
            if(${RUN_RESULT_VAR})
                set(${NAME} true)
            else()
                set(${NAME} false)
            endif()
        else()
            set(${NAME} ${RUN_RESULT_VAR})
        endif()
    endif()
endmacro()

macro(physica_cache_probe)
    physica_try_run(HwlocProbe HostDevAttr/HwlocProbe.cpp)
    if(NOT ${RUN_RESULT_VAR} EQUAL 0)
        message(FATAL_ERROR "Cache probe failed to run (${RUN_RESULT_VAR}):\n${COMPILE_LOG}\n${RUN_OUTPUT_VAR}")
    endif()
    string(REGEX MATCH "L1D=([0-9]+) L2=([0-9]+) L3=([0-9]+) line=([0-9]+)" _unused "${RUN_OUTPUT_VAR}")

    if("${CMAKE_MATCH_1}" STREQUAL "" OR "${CMAKE_MATCH_1}" EQUAL 0
       OR "${CMAKE_MATCH_2}" STREQUAL "" OR "${CMAKE_MATCH_2}" EQUAL 0
       OR "${CMAKE_MATCH_3}" STREQUAL "" OR "${CMAKE_MATCH_3}" EQUAL 0
       OR "${CMAKE_MATCH_4}" STREQUAL "" OR "${CMAKE_MATCH_4}" EQUAL 0)
        message(FATAL_ERROR "Cache probe output unparsable: ${RUN_OUTPUT_VAR}")
    endif()

    set(CacheSizeL1D ${CMAKE_MATCH_1} CACHE STRING "L1D data cache size in bytes")
    set(CacheSizeL2 ${CMAKE_MATCH_2} CACHE STRING "L2 cache size in bytes")
    set(CacheSizeL3 ${CMAKE_MATCH_3} CACHE STRING "L3 cache size in bytes")
    set(LineSizeL1D ${CMAKE_MATCH_4} CACHE STRING "L1D cache line size in bytes")
    mark_as_advanced(CacheSizeL1D CacheSizeL2 CacheSizeL3 LineSizeL1D)

    foreach(_var CacheSizeL1D CacheSizeL2 CacheSizeL3 LineSizeL1D)
        if("${${_var}}" STREQUAL "")
            message(FATAL_ERROR "${_var} is neither pinned nor probed: ${RUN_OUTPUT_VAR}")
        endif()

        if(${${_var}} EQUAL 0)
            message(FATAL_ERROR "Cache probe reported zero for ${_var}: ${RUN_OUTPUT_VAR}")
        endif()
    endforeach()
    message(STATUS "CacheProbe: L1D=${CacheSizeL1D} L2=${CacheSizeL2} L3=${CacheSizeL3} line=${LineSizeL1D}")
endmacro()
