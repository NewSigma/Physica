function(physica_add_test NAME SOURCES)
    add_executable(PhysicaTest-${NAME} ${SOURCES})
    add_test(NAME PhysicaTest-${NAME}
            COMMAND PhysicaTest-${NAME})
    target_link_libraries(PhysicaTest-${NAME} PhysicaCore)
endfunction()