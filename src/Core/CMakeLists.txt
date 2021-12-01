##############################################Libs################################################
#############################################Project##############################################
file(GLOB_RECURSE SOURCES_CORE *.cpp *.h *.cu)
add_library(PhysicaCore SHARED ${SOURCES_CORE})

install(TARGETS PhysicaCore DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/Core DESTINATION include/Physica)
#############################################Options##############################################
target_link_libraries(PhysicaCore PhysicaLogger PhysicaUtils pthread)
