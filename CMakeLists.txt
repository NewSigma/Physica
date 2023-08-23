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
option(PHYSICA_CUDA "Enable CUDA support" OFF)
option(PHYSICA_TORCH "Enable Torch support" OFF)
option(PHYSICA_HDF5 "Enable HDF5 support" OFF)
option(PHYSICA_OPTIMIZE "Enable aggressive optimizations" ON)
set(UseASM 1 CACHE STRING "Enable inline assembly code")
include(FindSIMD)

message(STATUS "Running under ${CMAKE_BUILD_TYPE} mode.")
# Policies
cmake_policy(SET CMP0074 NEW)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

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
    if (${CMAKE_BUILD_TYPE} MATCHES "RelWithDebugInfo")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -g -G -dopt")
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
-mrdrnd \
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
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
include_directories(3rdparty)

find_package(FFTW3 REQUIRED)
find_package(FFTW3f REQUIRED)
link_directories(${FFTW3_LIBRARY_DIRS})
link_directories(${FFTW3f_LIBRARY_DIRS})
include_directories(${FFTW3_INCLUDE_DIRS})

if(${PHYSICA_GUI})
    find_package(Qt6 COMPONENTS Core Widgets REQUIRED)
endif()

if(${PHYSICA_PROFILE})
    include_directories(${GPerfTools_INCLUDE_DIR})
    find_package(GPerfTools REQUIRED)
endif()

if(${PHYSICA_TORCH})
    find_package(Torch REQUIRED)
    include_directories(${TORCH_INCLUDE_DIRS})
endif()

if(${PHYSICA_HDF5})
    find_package(HDF5 NAMES hdf5 REQUIRED COMPONENTS CXX shared)
    include_directories(${HDF5_INCLUDE_DIR})
    add_definitions(-DPHYSICA_HDF5 -DH5_NO_DEPRECATED_SYMBOLS)
    string(REPLACE include lib HDF5_LIBRARY_DIR ${HDF5_INCLUDE_DIR}) # Cmake file of HDF5-1.14.1 seems do not provide HDF5_LIBRARY_DIR
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

if(${BUILD_BENCHMARKS})
    add_subdirectory(benchmark)
endif()

if(${PHYSICA_GUI})
    set(CMAKE_AUTOMOC ON) #target Physica needs MOC also, we put it after other targets to avoid unnecessary targets.(e.g. xxx_autogen)
    add_executable(Physica Physica.cpp ${HEADERS})
    target_link_libraries(Physica ${GPerfTools_LIBRARY})

    install(TARGETS Physica DESTINATION bin)
    #############################################Options##############################################
    target_link_libraries(Physica PhysicaCore PhysicaGui Qt6::Core Qt6::Widgets)
    if(${PHYSICA_TORCH})
        target_link_libraries(Physica PhysicaAI ${TORCH_LIBRARIES})
    endif()
    if(${PHYSICA_HDF5})
        target_link_libraries(Physica ${HDF5_CXX_SHARED_LIBRARY})
    endif()
endif()
install(DIRECTORY ${CMAKE_SOURCE_DIR}/cmake DESTINATION share)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/3rdparty/vectorclass DESTINATION include)
install(FILES include/Physica/Config.h DESTINATION include/Physica)
install(FILES include/Physica/SystemBits.h DESTINATION include/Physica)
