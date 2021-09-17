##############################################Libs################################################
#############################################Project##############################################
file(GLOB_RECURSE SOURCES_AI *.cpp *.h)
add_library(PhysicaAI SHARED ${SOURCES_AI})

install(TARGETS PhysicaAI DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/AI DESTINATION include/Physica)
#############################################Options##############################################
target_link_libraries(PhysicaAI PhysicaCore)
