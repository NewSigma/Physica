include(./PhysicaCheck.cmake)

find_package(Git QUIET)
if (${Git_FOUND})
    execute_process(COMMAND ${GIT_EXECUTABLE} log -1 --pretty=format:%H
                    OUTPUT_VARIABLE GIT_HASH
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

#CheckRightShift
physica_try_run(CheckRightShift Check/CheckRightShift.cpp)
if(${RUN_RESULT_VAR} MATCHES 1)
    message(STATUS "Performing CheckRightShift --Failed")
else()
    message(STATUS "Performing CheckRightShift --Passed")
endif()

physica_try_run(L1D HostDevAttr/L1D.cpp)
set(CacheSizeL1D ${RUN_OUTPUT_VAR})

physica_try_run(CacheLine HostDevAttr/CacheLine.cpp)
set(LineSizeL1D ${RUN_OUTPUT_VAR})

physica_try_run(L2 HostDevAttr/L2.cpp)
set(CacheSizeL2 ${RUN_OUTPUT_VAR})

physica_try_run(L3 HostDevAttr/L3.cpp)
set(CacheSizeL3 ${RUN_OUTPUT_VAR})

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
