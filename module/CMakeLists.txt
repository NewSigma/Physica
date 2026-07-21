add_subdirectory(Core)
add_subdirectory(Gui)
add_subdirectory(Logger)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/module DESTINATION module)
