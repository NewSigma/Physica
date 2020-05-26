cmake_minimum_required(VERSION 3.15)

project(Physica VERSION 0.0.1 LANGUAGES CXX CUDA)
#############################################Settings#############################################
option(RUN_CHECKS "Run Checks" OFF)
set(DEBUG_MODE ${CMAKE_BUILD_TYPE} MATCHES "Debug")
message(STATUS "Running under ${CMAKE_BUILD_TYPE} mode.")
#Platform related
set(GPU_ARCHITECTURE compute_75)
set(GPU_CODE sm_75)

if(${DEBUG_MODE})
    add_compile_options(-O0)
else()
    #FIXME: -On will cause segmentation fault.
    add_compile_options(-O2)
endif()
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)

set(CMAKE_CUDA_FLAGS "\
--compiler-options -fPIC \
--gpu-architecture ${GPU_ARCHITECTURE} \
--gpu-code ${GPU_CODE} \
--relocatable-device-code=true")

if(${RUN_CHECKS})
    add_subdirectory(Config)
endif()
##############################################Libs################################################
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
#############################################Project##############################################
include_directories(./)
add_subdirectory(AI)
add_subdirectory(Core)
add_subdirectory(Gui)
add_subdirectory(Interpreter)
if(${DEBUG_MODE})
    add_subdirectory(Test)
endif()

add_executable(Physica Physica.cpp)
#############################################Options##############################################
target_link_libraries(Physica Physica_Core Physica_Gui)
if(${DEBUG_MODE})
    target_link_libraries(Physica Physica_Test)
endif()