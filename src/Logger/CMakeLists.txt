#############################################Project##############################################
file(GLOB_RECURSE SOURCES_LOGGER *.cpp)
add_library(PhysicaLogger SHARED ${SOURCES_LOGGER})
#############################################Options##############################################