##############################################Libs################################################
#############################################Project##############################################
file(GLOB_RECURSE SOURCES_CORE *.cpp)
add_library(PhysicaUtils SHARED ${SOURCES_CORE})

install(TARGETS PhysicaUtils DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/Utils DESTINATION include/Physica)
#############################################Options##############################################
target_link_libraries(PhysicaUtils)
