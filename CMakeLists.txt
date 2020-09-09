cmake_minimum_required(VERSION 3.15)

project(Physica VERSION 0.0.1 LANGUAGES CXX CUDA)
#include all macros
include(cmake/Modules/PhysicaCheck.cmake)
#############################################Settings#############################################
option(RUN_CHECKS "Run Checks" OFF)
set(UseASM 1)
set(DEBUG_MODE ${CMAKE_BUILD_TYPE} MATCHES "Debug")
message(STATUS "Running under ${CMAKE_BUILD_TYPE} mode.")
#Platform related
set(GPU_ARCHITECTURE compute_75)
set(GPU_CODE sm_75)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED TRUE)
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)

if(${DEBUG_MODE})
    set(CMAKE_CXX_FLAGS -O0)
else()
    set(CMAKE_CXX_FLAGS -O2)
endif()

set(CMAKE_CXX_FLAGS "\
${CMAKE_CXX_FLAGS} \
-fPIC")

set(CMAKE_CUDA_FLAGS "\
--gpu-architecture ${GPU_ARCHITECTURE} \
--gpu-code ${GPU_CODE} \
--relocatable-device-code=true")

if(${RUN_CHECKS})
    add_subdirectory(config)
endif()
##############################################Libs################################################
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

find_package(Qt5 COMPONENTS Core Widgets REQUIRED)
#############################################Project##############################################
include_directories(./include)
file(GLOB_RECURSE HEADERS include/*.h)
add_subdirectory(src)
#Enable Tests
enable_testing()
add_subdirectory(test)

add_executable(Physica Physica.cpp ${HEADERS})
#############################################Options##############################################
target_link_libraries(Physica Physica_Core Physica_Gui Qt5::Core Qt5::Widgets)