##############################################Libs################################################
find_package(Qt5 COMPONENTS Core REQUIRED)
#############################################Project##############################################
file(GLOB_RECURSE SOURCES_INTERPRETER *.cpp *.h)
include_directories(Header)
add_library(Physica_Interpreter SHARED ${SOURCES_INTERPRETER})
#############################################Options##############################################
target_link_libraries(Physica_Interpreter Qt5::Core)
#Qt::Core will pass wrong option -fPIC to nvcc, so we have to remove it.
if (TARGET Qt5::Core)
    get_property(core_options TARGET Qt5::Core PROPERTY INTERFACE_COMPILE_OPTIONS)
    string(REPLACE "-fPIC" "" new_core_options "${core_options}")
    set_property(TARGET Qt5::Core PROPERTY INTERFACE_COMPILE_OPTIONS ${new_core_options})
endif()
target_compile_options(Physica_Interpreter PUBLIC $<$<COMPILE_LANGUAGE:CXX>:
        -O${OPTIMIZE_LEVEL}
        >)