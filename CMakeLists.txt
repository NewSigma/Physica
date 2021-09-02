cmake_minimum_required(VERSION 3.19)

project(Physica VERSION 0.0.1 LANGUAGES CXX CUDA)
#include all macros
include(cmake/Modules/PhysicaCheck.cmake)
#############################################Settings#############################################
option(BUILD_EXAMPLES "Build examples" ON)
option(BUILD_TESTS "Build tests" ON)
option(RUN_CHECKS "Run Checks" ON)
option(PHYSICA_PROFILE "Build with profiling support" OFF)
set(UseASM 1)
set(DEBUG_MODE ${CMAKE_BUILD_TYPE} MATCHES "Debug")
message(STATUS "Running under ${CMAKE_BUILD_TYPE} mode.")
# Policies
cmake_policy(SET CMP0074 NEW)
# Platform related
set(GPU_ARCHITECTURE compute_75)
set(GPU_CODE sm_75)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED TRUE)
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
set(CMAKE_CUDA_ARCHITECTURES 75)

if(${DEBUG_MODE})
    set(CMAKE_CXX_FLAGS -O0)
else()
    set(CMAKE_CXX_FLAGS -O2)
endif()

set(CMAKE_CXX_FLAGS "\
${CMAKE_CXX_FLAGS} \
-Wall -Wextra \
-fPIC \
-D__FILENAME__='\"$(basename $(notdir $(abspath $<)))\"'") #Define macro __FILENAME__

set(CMAKE_CUDA_FLAGS "\
${CMAKE_CUDA_FLAGS} \
--gpu-architecture ${GPU_ARCHITECTURE} \
--gpu-code ${GPU_CODE} \
--relocatable-device-code=true \
-D__FILENAME__='\"$(basename $(notdir $(abspath $<)))\"'")

if(${RUN_CHECKS})
    add_subdirectory(config)
endif()
##############################################Libs################################################
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
include_directories(3rdparty)

find_package(Qt5 COMPONENTS Core Widgets REQUIRED)
if(${PHYSICA_PROFILE})
    include_directories(${GPerfTools_INCLUDE_DIR})
    find_package(GPerfTools REQUIRED)
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

set(CMAKE_AUTOMOC ON) #target Physica needs MOC also, we put it after other targets to avoid unnecessary targets.(e.g. xxx_autogen)
add_executable(Physica Physica.cpp ${HEADERS})
target_link_libraries(Physica ${GPerfTools_LIBRARY})
target_compile_options(Physica PRIVATE -Wl,--rpath=${CMAKE_INSTALL_PREFIX}/lib)

install(TARGETS Physica DESTINATION bin)
install(FILES include/Physica/Config.h DESTINATION include/Physica)
install(FILES include/Physica/PhysicaInit.h DESTINATION include/Physica)
install(FILES include/Physica/SystemBits.h DESTINATION include/Physica)
#############################################Options##############################################
target_link_libraries(Physica PhysicaCore PhysicaGui Qt5::Core Qt5::Widgets)
