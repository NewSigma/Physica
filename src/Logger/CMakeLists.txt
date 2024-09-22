add_library(PhysicaLogger)

file(GLOB_RECURSE SOURCES_LOGGER *.cpp)
target_sources(PhysicaLogger PRIVATE ${SOURCES_LOGGER})

install(TARGETS PhysicaLogger DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/Logger DESTINATION include/Physica)
