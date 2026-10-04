add_subdirectory(Core)
if(${PHYSICA_GUI})
    add_subdirectory(Gui)
endif()
add_subdirectory(Logger)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/module DESTINATION module)
