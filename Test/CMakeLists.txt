#############################################Project##############################################
file(GLOB_RECURSE SOURCES_TEST *.cpp *.cu *.h *.cuh)
add_library(Physica_Test SHARED ${SOURCES_TEST})
#############################################Options##############################################
#Qt::Core will pass wrong option -fPIC to nvcc, so we have to remove it.
if ( TARGET Qt5::Core )
    get_property(core_options TARGET Qt5::Core PROPERTY INTERFACE_COMPILE_OPTIONS)
    string(REPLACE "-fPIC" "" new_core_options "${core_options}")
    set_property(TARGET Qt5::Core PROPERTY INTERFACE_COMPILE_OPTIONS ${new_core_options})
endif()
target_compile_options(Physica_Test PUBLIC $<$<COMPILE_LANGUAGE:CXX>:
        -O${OPTIMIZE_LEVEL}
        >)
target_compile_options(Physica_Test PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:
        --compiler-options -fPIC
        --gpu-architecture ${GPU_ARCHITECTURE}
        --gpu-code ${GPU_CODE}
        --optimize ${OPTIMIZE_LEVEL}
        --relocatable-device-code=true
        >)

target_link_libraries(Physica_Test Physica_AI Physica_Core)