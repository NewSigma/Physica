file(GLOB_RECURSE HEADERS ${CMAKE_HOME_DIRECTORY}/include/*.h ${CMAKE_HOME_DIRECTORY}/include/*.cuh)
set(PHYSICA_EXAMPLE_LIBS PhysicaCore pthread ${GPerfTools_LIBRARY})

if(${PHYSICA_TORCH})
set(PHYSICA_EXAMPLE_LIBS ${PHYSICA_EXAMPLE_LIBS} PhysicaAI ${TORCH_LIBRARIES})
endif()

if(${PHYSICA_GUI})
    set(CMAKE_AUTOMOC ON)
    set(PHYSICA_EXAMPLE_LIBS ${PHYSICA_EXAMPLE_LIBS} PhysicaGui)
endif()

function(physica_add_example NAME SOURCES)
    add_executable(Example-${NAME} ${SOURCES} ${HEADERS})
    target_link_libraries(Example-${NAME} ${PHYSICA_EXAMPLE_LIBS})
    install(TARGETS Example-${NAME} DESTINATION examples)
endfunction()
