cmake_minimum_required(VERSION 3.15)

project(Physica)
#############################################Settings#############################################
set(OPTIMIZE_LEVEL 0)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
##############################################Libs################################################
#############################################Project##############################################
add_subdirectory(Core)
add_subdirectory(Gui)
add_executable(Physica Physica.cpp Physica.h)
#############################################Options##############################################
target_compile_options(Physica PUBLIC $<$<COMPILE_LANGUAGE:CXX>:
        -O${OPTIMIZE_LEVEL}
        >)
target_link_libraries(Physica Physica_Core)
target_link_libraries(Physica Physica_Gui)