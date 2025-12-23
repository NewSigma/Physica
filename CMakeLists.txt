cmake_minimum_required(VERSION 3.24)

project(Physica
        VERSION 2026
        DESCRIPTION "Modern template library for computational physics"
        HOMEPAGE_URL https://gitee.com/newsigma/Physica
        LANGUAGES CXX)

if(${CMAKE_SOURCE_DIR} MATCHES ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In source build is a bad practice")
endif()

if(DEFINED SKBUILD_PROJECT_VERSION)
    if(NOT ${SKBUILD_PROJECT_VERSION} MATCHES ${CMAKE_PROJECT_VERSION})
        message(FATAL_ERROR "In consistent versions")
    endif()
endif()
message(STATUS "Running under ${CMAKE_BUILD_TYPE} mode.")
#############################################Settings#############################################
option(BUILD_TESTS "Build tests" ON)
option(BUILD_EXAMPLES "Build examples" OFF)
option(BUILD_BENCHMARKS "Build benchmarks" OFF)
option(BUILD_DOC "Build documentation" OFF)
option(RUN_CHECKS "Run Checks" ON)
option(PHYSICA_GUI "Enable GUI support" OFF)
option(PHYSICA_MKL "Enable MKL support" OFF)
option(PHYSICA_MPI "Enable MPI support" OFF)
option(PHYSICA_HDF5 "HDF5 is the standard IO format in Physica" ON)
option(PHYSICA_CUDA "Enables CUDA accleration" OFF)
option(PHYSICA_PYTHON "Enable python suppert" OFF)
option(PHYSICA_MIMALLOC "Use mimalloc" OFF)
option(PHYSICA_PROFILE "Build with profiling support" OFF)
option(PHYSICA_CCACHE ON)
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
# Dev options
option(PHYSICA_LLVMIR "Collect LLVM IR for benchmark" OFF)
option(PHYSICA_SANITIZE "Check memory error" OFF)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)
include(PhysicaConfigImpl)

if(${RUN_CHECKS})
    add_subdirectory(cmake/config)
endif()
add_subdirectory(cmake/Modules)
##############################################Project################################################
if (${PHYSICA_CCACHE} AND NOT ${PHYSICA_LLVMIR})
    find_program(CCACHE_PROGRAM ccache)
    if(CCACHE_PROGRAM)
        set(CMAKE_C_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
        set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
    else()
        message(STATUS "ccache not found")
    endif()
endif()

set(Physica_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")
include_directories(SYSTEM 3rdparty)
include_directories(${Physica_INCLUDE_DIR})
add_subdirectory(src)

if(${BUILD_TESTS})
    add_subdirectory(test)
endif()

if(${BUILD_EXAMPLES})
    if(NOT ${PHYSICA_GUI})
        message(FATAL_ERROR "To use Physica examples, you should enable PHYSICA_GUI.")
    endif()
    add_subdirectory(examples)
endif()

if(${BUILD_BENCHMARKS})
    if(NOT ${CMAKE_BUILD_TYPE} MATCHES "Release")
        message(WARNING "Running benchmark under non-release mode is not what you usually want, ignoring")
    else()
        add_subdirectory(benchmark)
    endif()
endif()

if(${BUILD_DOC})
    if(NOT ${PHYSICA_GUI})
        message(FATAL_ERROR "To use Physica documentation, you should enable PHYSICA_GUI.")
    endif()
    find_package(Sphinx REQUIRED)
    add_custom_target(Sphinx
                      COMMAND ${SPHINX_EXEC} -b html ${CMAKE_SOURCE_DIR}/doc ${CMAKE_BINARY_DIR}/sphinx
                      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                      COMMENT "Generating documentation with Sphinx")
    add_subdirectory(doc)
endif()

install(DIRECTORY ${CMAKE_SOURCE_DIR}/3rdparty/vectorclass DESTINATION include)
install(FILES include/Physica/Config.h DESTINATION include/Physica)
install(FILES include/Physica/CRCoro.h DESTINATION include/Physica)
install(FILES include/Physica/CRTPBase.h DESTINATION include/Physica)
install(FILES include/Physica/Macro.h DESTINATION include/Physica)
install(FILES include/Physica/PlainStruct.h DESTINATION include/Physica)
