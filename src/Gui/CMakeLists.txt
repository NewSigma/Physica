add_library(PhysicaGui)

file(GLOB_RECURSE SOURCES_GUI *.cpp)
target_sources(PhysicaGui PRIVATE ${SOURCES_GUI})
target_precompile_headers(PhysicaGui REUSE_FROM PhysicaCore)
target_link_libraries(PhysicaGui PUBLIC Physica::Gui)

install(TARGETS PhysicaGui DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/Gui DESTINATION include/Physica)
