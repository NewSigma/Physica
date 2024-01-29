cmake_minimum_required(VERSION 3.19)

project(Physica VERSION 0.0.1 LANGUAGES C CXX)
#include all macros
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)
include(PhysicaCheck)
#############################################Settings#############################################
option(BUILD_EXAMPLES "Build examples" ON)
option(BUILD_TESTS "Build tests" ON)
option(BUILD_BENCHMARKS "Build benchmarks" OFF)
option(RUN_CHECKS "Run Checks" ON)
option(PHYSICA_PROFILE "Build with profiling support" OFF)
option(PHYSICA_GUI "Enable GUI support" OFF)
option(PHYSICA_MKL "Enable MKL support" OFF)
option(PHYSICA_CUDA "Enable CUDA support" OFF)
option(PHYSICA_PYTHON "Enable python suppert" OFF)
option(PHYSICA_OPTIMIZE "Enable aggressive optimizations" ON)
set(UseASM 1 CACHE STRING "Enable inline assembly code")

message(STATUS "Running under ${CMAKE_BUILD_TYPE} mode.")
# Policies
cmake_policy(SET CMP0074 NEW)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")
set(CMAKE_VISIBILITY_INLINES_HIDDEN TRUE)
if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -G")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
                         -O2 \
                         -fno-math-errno \
                         -fno-signaling-nans \
                         -fno-trapping-math")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O2 -DNDEBUG")
    if (${CMAKE_BUILD_TYPE} MATCHES "RelWithDebInfo")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -g -G -dopt=on")
    endif()
endif()

if(${PHYSICA_OPTIMIZE})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

set(CMAKE_CXX_FLAGS "\
${CMAKE_CXX_FLAGS} \
-Wall \
-Wextra \
-fPIC \
-mrdrnd")

if(${PHYSICA_CUDA})
    if(NOT DEFINED PHYSICA_CUDA_ARCHITECTURES)
        message(FATAL_ERROR "To use CUDA, you should define PHYSICA_CUDA_ARCHITECTURES.")
    endif()
    add_definitions(-DPHYSICA_CUDA)
    find_package(CUDAToolkit REQUIRED)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED TRUE)
    set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
    set(CMAKE_CUDA_ARCHITECTURES ${PHYSICA_CUDA_ARCHITECTURES})
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
    add_subdirectory(config)
endif()
##############################################Libs################################################
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# vectorclass
include_directories(3rdparty)
# FFTW3
find_package(FFTW3 REQUIRED)
find_package(FFTW3f REQUIRED)
link_directories(${FFTW3_LIBRARY_DIRS})
link_directories(${FFTW3f_LIBRARY_DIRS})
include_directories(${FFTW3_INCLUDE_DIRS})
# HDF5
find_package(HDF5 NAMES hdf5 REQUIRED COMPONENTS CXX shared)
include_directories(${HDF5_INCLUDE_DIR})
add_definitions(-DPHYSICA_HDF5 -DH5_NO_DEPRECATED_SYMBOLS)
string(REPLACE include lib HDF5_LIBRARY_DIR ${HDF5_INCLUDE_DIR}) # Cmake file of HDF5-1.14.1 seems do not provide HDF5_LIBRARY_DIR

if(${PHYSICA_GUI})
    find_package(Qt6 COMPONENTS Core Widgets REQUIRED)
endif()

if(${PHYSICA_PROFILE})
    find_package(GPerfTools REQUIRED)
    include_directories(${GPerfTools_INCLUDE_DIR})
endif()

if(${PHYSICA_MKL})
    set(MKL_ARCH intel64)
    set(MKL_LINK dynamic)
    set(MKL_THREADING sequential)
    set(MKL_INTERFACE_FULL intel_ilp64)
    find_package(MKL REQUIRED)
    add_definitions(-DPHYSICA_MKL)
    include_directories(${MKL_INCLUDE})
endif()

if(${PHYSICA_CUDA})
    include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
endif()

if(${PHYSICA_PYTHON})
    set(PYBIND11_FINDPYTHON ON)
    add_subdirectory(3rdparty/pybind11)
    include_directories(${pybind11_INCLUDE_DIRS})
    add_definitions(-DPHYSICA_PYTHON)
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

if(NOT ${CMAKE_BUILD_TYPE} MATCHES "Debug")
    if(${BUILD_BENCHMARKS})
        add_subdirectory(benchmark)
    endif()
endif()

if(${PHYSICA_GUI})
    set(CMAKE_AUTOMOC ON) #target Physica needs MOC also, we put it after other targets to avoid unnecessary targets.(e.g. xxx_autogen)
    add_executable(Physica Physica.cpp ${HEADERS})
    target_link_libraries(Physica ${GPerfTools_LIBRARY})

    install(TARGETS Physica DESTINATION bin)
    #############################################Options##############################################
    target_link_libraries(Physica PhysicaCore PhysicaGui Qt6::Core Qt6::Widgets ${HDF5_CXX_SHARED_LIBRARY})
endif()
install(DIRECTORY ${CMAKE_SOURCE_DIR}/cmake DESTINATION share)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/3rdparty/vectorclass DESTINATION include)
install(FILES include/Physica/Config.h DESTINATION include/Physica)
install(FILES include/Physica/Macro.h DESTINATION include/Physica)
