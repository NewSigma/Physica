#############################################Project##############################################
file(GLOB_RECURSE SOURCES_TEST *.cpp *.cu *.h *.cuh)
add_library(Physica_Test SHARED ${SOURCES_TEST})
#############################################Options##############################################
target_link_libraries(Physica_Test Physica_AI Physica_Core)