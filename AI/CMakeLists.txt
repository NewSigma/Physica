##############################################Libs################################################
find_package(Qt5 COMPONENTS Core REQUIRED)
#############################################Project##############################################
file(GLOB_RECURSE SOURCES_AI *.cpp *.h)
include_directories(Header)
add_library(Physica_AI SHARED ${SOURCES_AI})
#############################################Options##############################################
target_link_libraries(Physica_AI Qt5::Core)

target_compile_options(Physica_AI PUBLIC $<$<COMPILE_LANGUAGE:CXX>:
        -O${OPTIMIZE_LEVEL}
        >)