add_library(PhysicaGui)

add_subdirectory(Plot)

target_link_libraries(PhysicaGui PUBLIC Physica::Gui)

set_target_properties(PhysicaGui PROPERTIES INSTALL_RPATH \${ORIGIN})

install(TARGETS PhysicaGui DESTINATION lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/Physica/Gui DESTINATION include/Physica)
