include(./PhysicaCheck.cmake)

find_package(Git QUIET)
if (${Git_FOUND})
    execute_process(COMMAND ${GIT_EXECUTABLE} log -1 --pretty=format:%H
                    OUTPUT_VARIABLE GIT_HASH
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

#TestBigEndian
include(TestBigEndian)
TEST_BIG_ENDIAN(PhysicaEndian)
if(${PhysicaEndian})
    set(PhysicaEndian PhysicaBigEndian)
else()
    set(PhysicaEndian PhysicaLittleEndian)
endif()
#CheckIntegerCast
physica_check(CheckIntegerCast Check/CheckIntegerCast.cpp)
if(${RUN_RESULT_VAR} MATCHES 1)
    message(STATUS "Performing CheckIntegerCast --Failed")
endif()
message(STATUS "Performing CheckIntegerCast --Passed")
#CheckRightShift
physica_check(CheckRightShift Check/CheckRightShift.cpp)
if(${RUN_RESULT_VAR} MATCHES 1)
    message(STATUS "Performing CheckRightShift --Failed")
endif()
message(STATUS "Performing CheckRightShift --Passed")

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
