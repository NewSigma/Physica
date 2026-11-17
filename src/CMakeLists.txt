include(HeaderCheck)

set(PHYSICA_COMPONENTS Core Logger)
add_subdirectory(Core)
add_subdirectory(Logger)

if(PHYSICA_GUI)
    add_subdirectory(Gui)
    list(APPEND PHYSICA_COMPONENTS Gui)
endif()

if(PHYSICA_TRANSFORMS)
    add_subdirectory(Transforms)
endif()

add_executable(Physica Physica.cpp)
target_link_libraries(Physica PhysicaCore ${GPerfTools_LIBRARY})
if(PHYSICA_GUI)
    target_link_libraries(Physica PhysicaGui)
endif()

foreach(component IN LISTS PHYSICA_COMPONENTS)
    install(EXPORT Physica${component}Targets FILE Physica${component}Targets.cmake NAMESPACE Physica:: DESTINATION "${PHYSICA_CMAKE_INSTALL_DIR}")
endforeach()
