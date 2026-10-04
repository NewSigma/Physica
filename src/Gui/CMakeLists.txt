add_library(PhysicaGui)

add_subdirectory(Plot)

target_link_libraries(PhysicaGui PUBLIC Physica::Gui)

set_target_properties(PhysicaGui PROPERTIES
    INSTALL_RPATH \${ORIGIN}
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR})

install(TARGETS PhysicaGui DESTINATION lib)
