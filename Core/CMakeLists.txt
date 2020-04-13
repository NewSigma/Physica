##############################################Libs################################################
find_package(Boost REQUIRED)
find_package(OpenCV REQUIRED)
find_package(Qt5 COMPONENTS Core REQUIRED)
#############################################Project##############################################
file(GLOB_RECURSE SOURCES_CORE *.cpp *.h)
include_directories(Header)
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
add_library(Physica_Core SHARED ${SOURCES_CORE})
#############################################Options##############################################
target_link_libraries(Physica_Core Qt5::Core)
#Qt::Core will pass wrong option -fPIC to nvcc, so we have to remove it.
if ( TARGET Qt5::Core )
    get_property(core_options TARGET Qt5::Core PROPERTY INTERFACE_COMPILE_OPTIONS)
    string(REPLACE "-fPIC" "" new_core_options "${core_options}")
    set_property(TARGET Qt5::Core PROPERTY INTERFACE_COMPILE_OPTIONS ${new_core_options})
endif()
target_compile_options(Physica_Core PUBLIC $<$<COMPILE_LANGUAGE:CXX>:
        -O${OPTIMIZE_LEVEL}
        >)
target_compile_options(Physica_Core PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:
        --compiler-options -fPIC
        --gpu-architecture ${GPU_ARCHITECTURE}
        --gpu-code ${GPU_ARCHITECTURE}
        --optimize ${GPU_CODE}
        --relocatable-device-code=true
        >)