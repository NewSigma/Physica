add_subdirectory(Core)
add_subdirectory(Logger)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/module DESTINATION module)
