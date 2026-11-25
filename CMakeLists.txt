# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2020-2026 Weibo He
# SPDX-License-Identifier: MPL-2.0
cmake_minimum_required(VERSION 3.29)

project(Physica
        VERSION 2026
        DESCRIPTION "Modern template library for computational physics"
        HOMEPAGE_URL https://gitee.com/newsigma/Physica
        LANGUAGES C CXX)

if(CMAKE_SOURCE_DIR MATCHES "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "In source build is a bad practice")
endif()

include(GNUInstallDirs)

message(STATUS "Running under ${CMAKE_BUILD_TYPE} mode.")
#############################################Settings#############################################
option(BUILD_TESTS "Build tests" ON)
option(BUILD_EXAMPLES "Build examples" OFF)
option(BUILD_BENCHMARKS "Build benchmarks" OFF)
option(BUILD_DOC "Build documentation" OFF)
option(RUN_CHECKS "Run Checks" ON)
option(PHYSICA_GUI "Lightweight Wrapper on Qt for 2D/3D Plotting" OFF)
option(PHYSICA_MKL "Enable MKL support" OFF)
option(PHYSICA_MPI "Enable MPI support" OFF)
option(PHYSICA_HDF5 "HDF5 is the standard IO format in Physica" OFF)
option(PHYSICA_CUDA "Enables CUDA acceleration" OFF)
option(PHYSICA_LLVM "" OFF)
option(PHYSICA_MIMALLOC "Use mimalloc" OFF)
# Dev options
set(PHYSICA_ARCH "native" CACHE STRING "Architecture for -march")
option(PHYSICA_PROFILE "Build with profiling support" OFF)
option(PHYSICA_TRANSFORMS "Enable custom LLVM passes" OFF)
option(PHYSICA_EmitLLVM "Collect LLVM IR for benchmark" OFF)
option(PHYSICA_ASAN "Compile with address sanitizer" OFF)
option(PHYSICA_NSAN "Compile with numerical sanitizer" OFF)
option(PHYSICA_CCACHE "" OFF)
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)
if(CMAKE_BUILD_TYPE MATCHES Debug)
    set(PHYSICA_TRANSFORMS OFF CACHE BOOL "Enable custom LLVM passes" FORCE)
endif()
set(Physica_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")
set(Physica_GENERATED_INCLUDE_DIR "${CMAKE_BINARY_DIR}/include")
include(PhysicaConfigImpl)

if(RUN_CHECKS)
    add_subdirectory(cmake/config)
endif()
set(PHYSICA_CMAKE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/Physica")
##############################################Project################################################
if(PHYSICA_CCACHE)
    find_program(CCACHE_PROGRAM ccache)
    if(CCACHE_PROGRAM)
        message(STATUS "Found ccache")
        set(CMAKE_C_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
        set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
        if(PHYSICA_CUDA)
            set(CMAKE_CUDA_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
        endif()
    else()
        message(STATUS "ccache not found")
    endif()
endif()

include_directories(SYSTEM 3rdparty)
include_directories(${Physica_INCLUDE_DIR})
include_directories(${Physica_GENERATED_INCLUDE_DIR})
add_subdirectory(src)

if(BUILD_TESTS)
    add_subdirectory(test)
endif()

if(BUILD_EXAMPLES)
    if(NOT PHYSICA_GUI)
        message(FATAL_ERROR "To use Physica examples, you should enable PHYSICA_GUI.")
    endif()
    add_subdirectory(examples)
endif()

if(BUILD_BENCHMARKS)
    if(NOT CMAKE_BUILD_TYPE MATCHES "Release")
        message(WARNING "Running benchmark under non-release mode is not what you usually want, ignoring")
    else()
        add_subdirectory(benchmark)
    endif()
endif()

if(BUILD_DOC)
    if(NOT PHYSICA_GUI)
        message(FATAL_ERROR "To use Physica documentation, you should enable PHYSICA_GUI.")
    endif()
    find_program(SPHINX_EXEC sphinx-build REQUIRED)
    add_custom_target(Sphinx
                      COMMAND ${SPHINX_EXEC} -b html -W ${CMAKE_SOURCE_DIR}/doc ${CMAKE_BINARY_DIR}/sphinx
                      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                      COMMENT "Generating documentation with Sphinx")
    add_subdirectory(doc)
endif()

##############################################Package###############################################
include(CMakePackageConfigHelpers)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/3rdparty/vectorclass DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES
            ${Physica_GENERATED_INCLUDE_DIR}/Physica/Config.h
            include/Physica/CRCoro.h
            include/Physica/CRTP.h
            include/Physica/Macro.h
            include/Physica/PlainStruct.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Physica)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/PhysicaConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/PhysicaConfig.cmake"
    INSTALL_DESTINATION "${PHYSICA_CMAKE_INSTALL_DIR}"
    PATH_VARS CMAKE_INSTALL_INCLUDEDIR)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/PhysicaConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)
install(FILES
            "${CMAKE_CURRENT_BINARY_DIR}/PhysicaConfig.cmake"
            "${CMAKE_CURRENT_BINARY_DIR}/PhysicaConfigVersion.cmake"
        DESTINATION "${PHYSICA_CMAKE_INSTALL_DIR}")

file(GLOB PHYSICA_MODULE_FILES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/*.cmake")
install(FILES ${PHYSICA_MODULE_FILES} DESTINATION "${PHYSICA_CMAKE_INSTALL_DIR}/modules")

install(DIRECTORY LICENSES DESTINATION "${CMAKE_INSTALL_DATADIR}")
install(FILES LICENSE DESTINATION "${CMAKE_INSTALL_DATADIR}")
