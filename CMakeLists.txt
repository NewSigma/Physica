cmake_minimum_required(VERSION 3.15)

project(Physica VERSION 0.0.1 LANGUAGES CXX CUDA)
#############################################Settings#############################################
option(RUN_CHECKS "Run Checks" OFF)
set(PHYSICA_COMPILE_TESTS True)
set(DEBUG_MODE ${CMAKE_BUILD_TYPE} MATCHES "Debug")
message(STATUS "Running under ${CMAKE_BUILD_TYPE} mode.")
#Platform related
set(GPU_ARCHITECTURE compute_75)
set(GPU_CODE sm_75)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)

if(${DEBUG_MODE})
    set(CMAKE_CXX_FLAGS -O0)
else()
    set(CMAKE_CXX_FLAGS -O2)
endif()

set(CMAKE_CUDA_FLAGS "\
--compiler-options -fPIC \
--gpu-architecture ${GPU_ARCHITECTURE} \
--gpu-code ${GPU_CODE} \
--relocatable-device-code=true")

if(${RUN_CHECKS})
    add_subdirectory(config)
endif()
#############################################Project##############################################
if(${PHYSICA_COMPILE_TESTS})
    add_subdirectory(test)
endif()

include_directories(./include)
file(GLOB_RECURSE HEADERS include/*.h)
add_subdirectory(src)
add_executable(Physica Physica.cpp ${HEADERS})
#############################################Options##############################################
target_link_libraries(Physica Physica_Core Physica_Gui)
if(${PHYSICA_COMPILE_TESTS})
    target_link_libraries(Physica Physica_Test)
endif()