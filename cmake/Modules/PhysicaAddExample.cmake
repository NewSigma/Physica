set(CMAKE_AUTOMOC ON)

file(GLOB_RECURSE HEADERS ${CMAKE_HOME_DIRECTORY}/include/*.h ${CMAKE_HOME_DIRECTORY}/include/*.cuh)

function(physica_add_example NAME SOURCES)
    add_executable(Example-${NAME} ${SOURCES} ${HEADERS})
    target_link_libraries(Example-${NAME} PhysicaCore PhysicaGui)
endfunction()
