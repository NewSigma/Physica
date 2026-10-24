include(./PhysicaCheck.cmake)

find_package(Git QUIET REQUIRED)
execute_process(COMMAND ${GIT_EXECUTABLE} log -1 --pretty=format:%H
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                OUTPUT_VARIABLE GIT_HASH
                ERROR_VARIABLE GIT_HASH_ERROR
                RESULT_VARIABLE GIT_HASH_RESULT
                OUTPUT_STRIP_TRAILING_WHITESPACE)
if (NOT GIT_HASH_RESULT EQUAL 0)
    message(FATAL_ERROR "Failed to read the git revision (git returned ${GIT_HASH_RESULT}): ${GIT_HASH_ERROR}")
endif()
if (GIT_HASH STREQUAL "")
    message(FATAL_ERROR "Git revision is empty?")
endif()

#CheckRightShift
physica_try_run(CheckRightShift Check/CheckRightShift.cpp)
if(${RUN_RESULT_VAR} MATCHES 1)
    message(STATUS "Performing CheckRightShift --Failed")
else()
    message(STATUS "Performing CheckRightShift --Passed")
endif()

find_package(Hwloc REQUIRED)
physica_try_run(HwlocProbe HostDevAttr/HwlocProbe.cpp)
if(NOT ${RUN_RESULT_VAR} EQUAL 0)
    message(FATAL_ERROR "Cache probe (hwloc) failed to run (${RUN_RESULT_VAR}):\n${COMPILE_LOG}\n${RUN_OUTPUT_VAR}")
endif()

string(REGEX MATCH "L1D=([0-9]+) L2=([0-9]+) L3=([0-9]+) line=([0-9]+)" _unused "${RUN_OUTPUT_VAR}")
set(CacheSizeL1D ${CMAKE_MATCH_1})
set(CacheSizeL2 ${CMAKE_MATCH_2})
set(CacheSizeL3 ${CMAKE_MATCH_3})
set(LineSizeL1D ${CMAKE_MATCH_4})
foreach(_var CacheSizeL1D CacheSizeL2 CacheSizeL3 LineSizeL1D)
    if("${${_var}}" STREQUAL "")
        message(FATAL_ERROR "Cache probe (hwloc) output unparsable (${_var}): ${RUN_OUTPUT_VAR}")
    endif()

    if(${${_var}} EQUAL 0)
        message(FATAL_ERROR "Cache probe (hwloc) reported zero for ${_var}: ${RUN_OUTPUT_VAR}")
    endif()
endforeach()
message(STATUS "CacheProbe(hwloc): L1D=${CacheSizeL1D} L2=${CacheSizeL2} L3=${CacheSizeL3} line=${LineSizeL1D}")

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
configure_file(Config.h.in ${CMAKE_SOURCE_DIR}/include/Physica/Config.h)
