cmake_minimum_required(VERSION 3.10)

project(Physica VERSION 0.0.1)
#############################################Settings#############################################
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
##############################################Libs################################################
enable_language(CUDA)
find_package(OpenCV REQUIRED)
find_package(Boost REQUIRED)
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
#############################################Project##############################################
add_subdirectory(Gui)
file(GLOB_RECURSE SOURCES Physica.cpp Kernel.cu Core/*.* Utils/*.*)
file(GLOB_RECURSE RESOURCES Resources/*.*)
add_executable(${PROJECT_NAME} ${SOURCES} ${RESOURCES})
#############################################Options##############################################
set(OPTIMIZE_LEVEL 0)
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_Gui)
target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS})
target_compile_options(${PROJECT_NAME} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:
        -O${OPTIMIZE_LEVEL}
        >)
target_compile_options(${PROJECT_NAME} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
        --gpu-architecture=sm_50
        --optimize ${OPTIMIZE_LEVEL}
        --use_fast_math
        --relocatable-device-code=true
        >)