cmake_minimum_required(VERSION 3.19)

project(Physica
        VERSION 0.0.1
        DESCRIPTION "Modern template library for computational physics"
        HOMEPAGE_URL https://gitee.com/newsigma/Physica
        LANGUAGES C CXX)
#include all macros
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)
#############################################Settings#############################################
option(BUILD_EXAMPLES "Build examples" OFF)
option(BUILD_TESTS "Build tests" ON)
option(BUILD_BENCHMARKS "Build benchmarks" OFF)
option(BUILD_DOC "Build documentation" OFF)
option(RUN_CHECKS "Run Checks" ON)
option(PHYSICA_PROFILE "Build with profiling support" OFF)
option(PHYSICA_GUI "Enable GUI support" OFF)
option(PHYSICA_MKL "Enable MKL support" OFF)
option(PHYSICA_CUDA "Enable CUDA support" OFF)
option(PHYSICA_PYTHON "Enable python suppert" OFF)
option(PHYSICA_OPTIMIZE "Enable aggressive optimizations" ON)
set(UseASM 1 CACHE STRING "Enable inline assembly code")
# Policies
cmake_policy(SET CMP0074 NEW)

message(STATUS "Running under ${CMAKE_BUILD_TYPE} mode.")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_VISIBILITY_PRESET "hidden")
set(CMAKE_VISIBILITY_INLINES_HIDDEN TRUE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fPIC -mrdrnd")
if(${PHYSICA_OPTIMIZE})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()
if(NOT ${CMAKE_BUILD_TYPE} MATCHES "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-math-errno -fno-signaling-nans -fno-trapping-math")
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-ignored-optimization-argument -Wno-unknown-warning-option")
    add_link_options(-lstdc++) # Add this if you prefer libstdc++
endif()

if(${PHYSICA_CUDA})
    if(NOT DEFINED PHYSICA_CUDA_ARCHITECTURES)
        message(FATAL_ERROR "To use CUDA, you should define PHYSICA_CUDA_ARCHITECTURES.")
    endif()
    add_definitions(-DPHYSICA_CUDA)
    find_package(CUDAToolkit REQUIRED)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD ${CMAKE_CXX_STANDARD})
    set(CMAKE_CUDA_STANDARD_REQUIRED ${CMAKE_CXX_STANDARD_REQUIRED})
    set(CMAKE_CUDA_VISIBILITY_PRESET ${CMAKE_CXX_VISIBILITY_PRESET})
    set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
    set(CMAKE_CUDA_ARCHITECTURES ${PHYSICA_CUDA_ARCHITECTURES})

    set(CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} -G")
    set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "${CMAKE_CUDA_FLAGS_RELWITHDEBINFO} -G -dopt=on")
    # --diag-suppress 186: Comparation between unsigned int and 0 in InverseMatrix does not exist
    set(CMAKE_CUDA_FLAGS "\
        ${CMAKE_CUDA_FLAGS} \
        --Werror cross-execution-space-call \
        --Wreorder \
        --Wdefault-stream-launch \
        --Wext-lambda-captures-this \
        --Wno-deprecated-declarations \
        --Wno-deprecated-gpu-targets \
        --relocatable-device-code=true \
        --diag-suppress 186 \
        --diag-suppress 20011 \
        --default-stream per-thread \
        ${CMAKE_CXX_FLAGS}")
endif()

if(${RUN_CHECKS})
    add_subdirectory(cmake/config)
endif()
##############################################Libs################################################
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# vectorclass
include_directories(SYSTEM 3rdparty)
# FFTW3
find_package(FFTW3 REQUIRED)
find_package(FFTW3f REQUIRED)
link_directories(${FFTW3_LIBRARY_DIRS})
link_directories(${FFTW3f_LIBRARY_DIRS})
include_directories(SYSTEM ${FFTW3_INCLUDE_DIRS})
# HDF5
find_package(HDF5 NAMES hdf5 REQUIRED COMPONENTS CXX shared)
include_directories(SYSTEM ${HDF5_INCLUDE_DIR_CPP})
add_definitions(-DPHYSICA_HDF5 -DH5_NO_DEPRECATED_SYMBOLS)

if(${PHYSICA_GUI})
    find_package(Qt6 COMPONENTS Core Widgets REQUIRED)
endif()

if(${PHYSICA_PROFILE})
    find_package(GPerfTools REQUIRED)
    include_directories(SYSTEM ${GPerfTools_INCLUDE_DIR})
endif()

if(${PHYSICA_MKL})
    set(MKL_ARCH intel64)
    set(MKL_LINK dynamic)
    set(MKL_THREADING sequential)
    set(MKL_INTERFACE_FULL intel_ilp64)
    find_package(MKL REQUIRED)
    add_definitions(-DPHYSICA_MKL)
    include_directories(SYSTEM ${MKL_INCLUDE})
endif()

if(${PHYSICA_CUDA})
    include_directories(SYSTEM ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
endif()

if(${PHYSICA_PYTHON})
    set(PYBIND11_FINDPYTHON ON)
    add_subdirectory(3rdparty/pybind11)
    include_directories(${pybind11_INCLUDE_DIRS})
    add_definitions(-DPHYSICA_PYTHON)
endif()

if(${BUILD_DOC})
    find_package(Sphinx REQUIRED)
endif()
#############################################Project##############################################
include_directories(./include)
file(GLOB_RECURSE HEADERS include/*.h include/*.cuh)
add_subdirectory(src)
#Enable examples
if(${BUILD_EXAMPLES})
    if(NOT ${PHYSICA_GUI})
        message(FATAL_ERROR "To use Physica examples, you should enable PHYSICA_GUI.")
    endif()
    add_subdirectory(examples)
endif()
#Enable Tests
if(${BUILD_TESTS})
    add_subdirectory(test)
endif()

if(${BUILD_BENCHMARKS})
    if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
        message(WARNING "Running benchmark under debug mode is not what you usually want, ignoring")
    else()
        add_subdirectory(benchmark)
    endif()
endif()

if(${BUILD_DOC})
    add_custom_target(Sphinx
                      COMMAND ${SPHINX_EXEC} -b html ${CMAKE_SOURCE_DIR}/doc ${CMAKE_BINARY_DIR}/sphinx
                      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                      COMMENT "Generating documentation with Sphinx")
endif()

install(DIRECTORY ${CMAKE_SOURCE_DIR}/cmake/Modules DESTINATION share/cmake)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/3rdparty/vectorclass DESTINATION include)
install(FILES include/Physica/Config.h DESTINATION include/Physica)
install(FILES include/Physica/Macro.h DESTINATION include/Physica)
