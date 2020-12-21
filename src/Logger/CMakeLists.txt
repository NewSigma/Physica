#############################################Project##############################################
file(GLOB_RECURSE SOURCES_LOGGER *.cpp)
add_library(PhysicaLogger SHARED ${SOURCES_LOGGER})

install(TARGETS PhysicaLogger DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/Logger DESTINATION include/Physica)
#############################################Options##############################################