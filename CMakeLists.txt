cmake_minimum_required(VERSION 3.15)

project(Physica VERSION 0.0.1 LANGUAGES CXX CUDA)
#############################################Settings#############################################
set(OPTIMIZE_LEVEL 0)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(ENABLE_TEST True)
#Platform related
set(GPU_ARCHITECTURE compute_75)
set(GPU_CODE sm_75)
##############################################Libs################################################
#############################################Project##############################################
include_directories(./)
add_subdirectory(Core)
include_directories(Core/Header)
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
add_subdirectory(Gui)
include_directories(Gui/Header)
if(ENABLE_TEST)
    add_subdirectory(Test)
    include_directories(Test)
endif()
add_executable(Physica Physica.cpp)
#############################################Options##############################################
target_compile_options(Physica PUBLIC $<$<COMPILE_LANGUAGE:CXX>:
        -O${OPTIMIZE_LEVEL}
        >)
target_link_libraries(Physica Physica_Core Physica_Gui)