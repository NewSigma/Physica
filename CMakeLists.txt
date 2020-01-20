cmake_minimum_required(VERSION 3.10)

project(Physica VERSION 0.0.1)
#############################################Settings#############################################
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_COMPILER /usr/local/cuda-10.1/bin/nvcc)
set(CMAKE_PREFIX_PATH "~/Qt5.13.2/5.13.2/gcc_64/lib/cmake")
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)
##############################################Libs################################################
enable_language(CUDA)
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
find_package(Qt5 COMPONENTS Core Gui Widgets REQUIRED)
find_package(OpenCV REQUIRED)
#############################################Project##############################################
if (TARGET Qt5::Core)
    get_property(core_options TARGET Qt5::Core PROPERTY INTERFACE_COMPILE_OPTIONS)
    string(REPLACE "-fPIC" "" new_core_options "${core_options}")
    set_property(TARGET Qt5::Core PROPERTY INTERFACE_COMPILE_OPTIONS ${new_core_options})
    set_property(TARGET Qt5::Core PROPERTY INTERFACE_POSITION_INDEPENDENT_CODE "ON")
    set(CMAKE_CXX_COMPILE_OPTIONS_PIE "-fPIC")
endif()

file(GLOB_RECURSE SOURCES Core/*.* Gui/*.* Utils/*.*)
file(GLOB_RECURSE RESOURCES Resources/*.*)
add_executable(${PROJECT_NAME} Physica.cpp Kernel.cu ${SOURCES} ${RESOURCES})
#############################################Options##############################################
target_link_libraries(${PROJECT_NAME} Qt5::Core Qt5::Gui Qt5::Widgets)
target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS})
target_compile_options(${PROJECT_NAME} PRIVATE $<$<COMPILE_LANGUAGE:CXX>:
        -O0
        >)
target_compile_options(${PROJECT_NAME} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
        --gpu-architecture=sm_35
        --use_fast_math
        --relocatable-device-code=true
        >)