include(./PhysicaCheck.cmake)

find_package(Git QUIET)
if (${Git_FOUND})
    execute_process(COMMAND ${GIT_EXECUTABLE} log -1 --pretty=format:%H
                    OUTPUT_VARIABLE GIT_HASH
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

#CheckRightShift
physica_check(CheckRightShift Check/CheckRightShift.cpp)
if(${RUN_RESULT_VAR} MATCHES 1)
    message(STATUS "Performing CheckRightShift --Failed")
endif()
message(STATUS "Performing CheckRightShift --Passed")

physica_check(L1D HostDevAttr/L1D.cpp)
set(CacheSizeL1D ${RUN_OUTPUT_VAR})

physica_check(CacheLine HostDevAttr/CacheLine.cpp)
set(LineSizeL1D ${RUN_OUTPUT_VAR})

physica_check(L2 HostDevAttr/L2.cpp)
set(CacheSizeL2 ${RUN_OUTPUT_VAR})

physica_check(L3 HostDevAttr/L3.cpp)
set(CacheSizeL3 ${RUN_OUTPUT_VAR})

if(${PHYSICA_CUDA})
    # WarpSize
    physica_check(WarpSize CUDADevAttr/WarpSize.cpp)
    set(WarpSize ${RUN_RESULT_VAR})
    message(STATUS "Performing WarpSize --${RUN_RESULT_VAR}")
    # MemoryPoolsSupported
    physica_check(MemoryPoolsSupported CUDADevAttr/MemoryPoolsSupported.cpp)
    if(${RUN_RESULT_VAR} MATCHES 0)
        set(MemoryPoolsSupported false)
    else()
        set(MemoryPoolsSupported true)
    endif()
    message(STATUS "Performing MemoryPoolsSupported --${MemoryPoolsSupported}")
endif()

message(STATUS)
configure_file(Config.h.in ${CMAKE_SOURCE_DIR}/include/Physica/Config.h)
