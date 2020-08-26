function(physica_add_test NAME SOURCES)
    add_executable(PhysicaTest_${NAME} ${SOURCES})
    add_test(NAME PhysicaTest_${NAME}
            COMMAND PhysicaTest_${NAME})
    target_link_libraries(PhysicaTest_${NAME} Physica_Core)
endfunction()