include(PhysicaTargets)

add_subdirectory(Core)

if(${PHYSICA_GUI})
    add_subdirectory(Gui)
endif()
add_subdirectory(Logger)    
if(${PHYSICA_PYTHON})
    add_subdirectory(Python)
endif()

add_executable(Physica Physica.cpp)
target_link_libraries(Physica PhysicaCore ${GPerfTools_LIBRARY})
if(${PHYSICA_GUI})
    target_link_libraries(Physica PhysicaGui)
    set_target_properties(Physica PROPERTIES AUTOMOC ON)
endif()
