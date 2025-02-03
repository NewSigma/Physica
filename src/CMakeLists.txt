include(PhysicaTargets)

add_subdirectory(Core)

if(${PHYSICA_GUI})
    add_subdirectory(Gui)
endif()
add_subdirectory(Logger)    
if(${PHYSICA_PYTHON})
    add_subdirectory(Python)
endif()

if(NOT DEFINED SKBUILD_PROJECT_VERSION)
    add_executable(Physica Physica.cpp)
    target_link_libraries(Physica PhysicaCore ${GPerfTools_LIBRARY})
    if(${PHYSICA_GUI})
        target_link_libraries(Physica PhysicaGui)
    endif()
endif()
