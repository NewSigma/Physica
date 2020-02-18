##############################################Libs################################################
find_package(CUDA)
find_package(Boost REQUIRED)
find_package(OpenCV REQUIRED)
#############################################Project##############################################
file(GLOB_RECURSE SOURCES_CORE *.*)
add_library(Physica_Core SHARED ${SOURCES_CORE})
#############################################Options##############################################
target_compile_options(Physica_Core PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
        --gpu-architecture=sm_50
        --optimize ${OPTIMIZE_LEVEL}
        --use_fast_math
        --relocatable-device-code=true
        >)