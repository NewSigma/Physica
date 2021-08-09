function(physica_add_test NAME SOURCES)
    add_executable(Test-${NAME} ${SOURCES})
    add_test(NAME Test-${NAME}
            COMMAND Test-${NAME})
    target_link_libraries(Test-${NAME} PhysicaCore ${GPerfTools_LIBRARY})
endfunction()