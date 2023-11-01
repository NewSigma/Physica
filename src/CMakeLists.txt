#############################################Project##############################################
add_subdirectory(Core)
if(${PHYSICA_GUI})
    add_subdirectory(Gui)
endif()
add_subdirectory(Logger)
add_subdirectory(Utils)