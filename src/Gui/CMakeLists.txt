add_library(PhysicaGui)

file(GLOB_RECURSE SOURCES_GUI *.cpp)
target_sources(PhysicaGui PRIVATE ${SOURCES_GUI})
target_link_libraries(PhysicaGui PUBLIC Physica::Gui)

set_target_properties(PhysicaGui PROPERTIES INSTALL_RPATH \${ORIGIN})

install(TARGETS PhysicaGui DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/Gui DESTINATION include/Physica)
