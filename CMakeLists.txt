cmake_minimum_required(VERSION 3.19)

project(Physica VERSION 0.0.1 LANGUAGES CXX)
#include all macros
include(cmake/Modules/PhysicaCheck.cmake)
#############################################Settings#############################################
option(BUILD_EXAMPLES "Build examples" ON)
option(BUILD_TESTS "Build tests" ON)
option(RUN_CHECKS "Run Checks" ON)
option(PHYSICA_PROFILE "Build with profiling support" OFF)
option(PHYSICA_FFTW3 "Use fftw3" OFF)
option(PHYSICA_GUI "Enable GUI support" OFF)
option(PHYSICA_CUDA "Enable CUDA support" OFF)
option(PHYSICA_TORCH "Enable Torch support" OFF)
set(PHYSICA_SIMD_OPTION "" CACHE STRING "Compiler options of enabled SIMD instruction sets") # Pass options based on processor to gain optimized performance
set(UseASM 1 CACHE STRING "Enable inline assembly code")

set(DEBUG_MODE ${CMAKE_BUILD_TYPE} MATCHES "Debug")
message(STATUS "Running under ${CMAKE_BUILD_TYPE} mode.")
# Policies
cmake_policy(SET CMP0074 NEW)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

if(${DEBUG_MODE})
    set(CMAKE_CXX_FLAGS -O0)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -G")
else()
    set(CMAKE_CXX_FLAGS -O2)
endif()

set(CMAKE_CXX_FLAGS "\
${CMAKE_CXX_FLAGS} \
${PHYSICA_SIMD_OPTION} \
-Wall -Wextra \
-fPIC \
-D__FILENAME__='\"$(basename $(notdir $(abspath $<)))\"'") #Define macro __FILENAME__

if(${PHYSICA_CUDA})
    add_definitions(-DPHYSICA_CUDA)
    find_package(CUDAToolkit REQUIRED)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED TRUE)
    set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
    set(CMAKE_CUDA_ARCHITECTURES 75)
    set(CMAKE_CUDA_FLAGS "\
    ${CMAKE_CUDA_FLAGS} \
    --relocatable-device-code=true \
    --diag-suppress 20011 \
    -D__FILENAME__='\"$(basename $(notdir $(abspath $<)))\"'")
endif()

if(${RUN_CHECKS})
    add_subdirectory(config)
endif()
##############################################Libs################################################
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
include_directories(3rdparty)

if(${PHYSICA_FFTW3})
    add_compile_definitions(PHYSICA_FFTW3)
    find_package(FFTW3 COMPONENTS fftw3 REQUIRED)
endif()

if(${PHYSICA_GUI})
    find_package(Qt6 COMPONENTS Core Widgets REQUIRED)
endif()

if(${PHYSICA_PROFILE})
    include_directories(${GPerfTools_INCLUDE_DIR})
    find_package(GPerfTools REQUIRED)
endif()

if(${PHYSICA_TORCH})
    find_package(Torch REQUIRED)
    include_directories(${Torch_INCLUDE_DIR})
endif()
#############################################Project##############################################
include_directories(./include)
file(GLOB_RECURSE HEADERS include/*.h include/*.cuh)
add_subdirectory(src)
#Enable examples
if(${BUILD_EXAMPLES})
    add_subdirectory(examples)
endif()
#Enable Tests
if(${BUILD_TESTS})
    add_subdirectory(test)
endif()

if(${PHYSICA_GUI})
    set(CMAKE_AUTOMOC ON) #target Physica needs MOC also, we put it after other targets to avoid unnecessary targets.(e.g. xxx_autogen)
    add_executable(Physica Physica.cpp ${HEADERS})
    target_link_libraries(Physica ${GPerfTools_LIBRARY})
    target_compile_options(Physica PRIVATE -Wl,--rpath=${CMAKE_INSTALL_PREFIX}/lib)

    install(TARGETS Physica DESTINATION bin)
    install(FILES include/Physica/Config.h DESTINATION include/Physica)
    install(FILES include/Physica/PhysicaInit.h DESTINATION include/Physica)
    install(FILES include/Physica/SystemBits.h DESTINATION include/Physica)
    #############################################Options##############################################
    target_link_libraries(Physica PhysicaCore PhysicaGui Qt6::Core Qt6::Widgets)
endif()
