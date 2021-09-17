##############################################Libs################################################
#############################################Project##############################################
file(GLOB_RECURSE SOURCES_INTERPRETER *.cpp *.h)
add_library(PhysicaInterpreter SHARED ${SOURCES_INTERPRETER})

install(TARGETS PhysicaInterpreter DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/Interpreter DESTINATION include/Physica)
#############################################Options##############################################
target_link_libraries(PhysicaInterpreter PhysicaCore)