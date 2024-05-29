set(PHYSICA_TEST_LIBS PhysicaCore ${GPerfTools_LIBRARY})

function(physica_add_test NAME SOURCES)
    add_executable(Test-${NAME} ${SOURCES} ${ARGN})
    add_test(NAME Test-${NAME} COMMAND Test-${NAME})
    target_link_libraries(Test-${NAME} ${PHYSICA_TEST_LIBS})
endfunction()
