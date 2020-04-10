#############################################Project##############################################
file(GLOB_RECURSE SOURCES_TEST Test.cpp *.h)
add_library(Physica_Test SHARED ${SOURCES_TEST})
#############################################Options##############################################
target_link_libraries(Physica_Test Physica_Core)