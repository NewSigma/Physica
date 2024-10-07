add_library(PhysicaGui)

set_target_properties(PhysicaGui PROPERTIES AUTOMOC ON)
set_target_properties(PhysicaGui PROPERTIES AUTOUIC ON)
set_target_properties(PhysicaGui PROPERTIES AUTORCC ON)

file(GLOB_RECURSE SOURCES_GUI *.cpp)
target_sources(PhysicaGui PRIVATE ${SOURCES_GUI})
target_precompile_headers(PhysicaGui REUSE_FROM PhysicaCore)
target_link_libraries(PhysicaGui PUBLIC Physica::Gui)

install(TARGETS PhysicaGui DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/Gui DESTINATION include/Physica)
