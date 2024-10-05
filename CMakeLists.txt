cmake_minimum_required(VERSION 3.19)

project(Physica
        VERSION 0.0.1
        DESCRIPTION "Modern template library for computational physics"
        HOMEPAGE_URL https://gitee.com/newsigma/Physica
        LANGUAGES C CXX)

if (${CMAKE_SOURCE_DIR} MATCHES ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In source build is a bad practice")
endif()
#############################################Settings#############################################
option(BUILD_TESTS "Build tests" ON)
option(BUILD_EXAMPLES "Build examples" OFF)
option(BUILD_BENCHMARKS "Build benchmarks" OFF)
option(BUILD_DOC "Build documentation" OFF)
option(RUN_CHECKS "Run Checks" ON)
option(PHYSICA_PROFILE "Build with profiling support" OFF)
option(PHYSICA_GUI "Enable GUI support" OFF)
option(PHYSICA_MKL "Enable MKL support" OFF)
option(PHYSICA_MPI "Enable MPI support" OFF)
option(PHYSICA_HDF5 "HDF5 is the standard IO format in Physica" ON)
option(PHYSICA_CUDA "CUDA enables GPU accleration" OFF)
option(PHYSICA_PYTHON "Enable python suppert" OFF)
option(PHYSICA_OPTIMIZE "Enable aggressive optimizations" ON)
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(CMAKE_DISABLE_PRECOMPILE_HEADERS "Set OFF to acclerate compile using PCH" ON)

message(STATUS "Running under ${CMAKE_BUILD_TYPE} mode.")
# Policies
cmake_policy(SET CMP0074 NEW)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)
include(PhysicaConfigImpl)

if(${RUN_CHECKS})
    add_subdirectory(cmake/config)
endif()
add_subdirectory(cmake/Modules)
##############################################Project################################################
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(Physica_INCLUDE_DIR "./include")
include_directories(SYSTEM 3rdparty)
include_directories(${Physica_INCLUDE_DIR})
file(GLOB_RECURSE HEADERS include/*.h include/*.cuh)
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
    if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
        message(WARNING "Running benchmark under debug mode is not what you usually want, ignoring")
    else()
        add_subdirectory(benchmark)
    endif()
endif()

if(${BUILD_DOC})
    find_package(Sphinx REQUIRED)
    add_custom_target(Sphinx
                      COMMAND ${SPHINX_EXEC} -b html ${CMAKE_SOURCE_DIR}/doc ${CMAKE_BINARY_DIR}/sphinx
                      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                      COMMENT "Generating documentation with Sphinx")
endif()

install(DIRECTORY ${CMAKE_SOURCE_DIR}/3rdparty/vectorclass DESTINATION include)
install(FILES include/Physica/Config.h DESTINATION include/Physica)
install(FILES include/Physica/CRTPBase.h DESTINATION include/Physica)
install(FILES include/Physica/Macro.h DESTINATION include/Physica)
install(FILES include/Physica/PlainStruct.h DESTINATION include/Physica)
