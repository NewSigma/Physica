function(physica_add_example NAME SOURCES)
    add_executable(Example-${NAME})
    target_sources(Example-${NAME} PRIVATE ${SOURCES} ${ARGN})
    target_link_libraries(Example-${NAME} PhysicaCore ${GPerfTools_LIBRARY})

    if(${PHYSICA_GUI})
        target_link_libraries(Example-${NAME} PhysicaGui)
    endif()
    install(TARGETS Example-${NAME} DESTINATION examples)
endfunction()

add_subdirectory(Core)
add_subdirectory(Gui)
