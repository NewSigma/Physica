configure_file(PhysicaConfig.cmake.in ${CMAKE_BINARY_DIR}/PhysicaConfig.cmake @ONLY)

install(FILES ${CMAKE_BINARY_DIR}/PhysicaConfig.cmake DESTINATION share/cmake/Modules)
file(GLOB ModuleFiles *.cmake)
install(FILES ${ModuleFiles} DESTINATION share/cmake/Modules)
