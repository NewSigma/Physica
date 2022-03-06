function(physica_add_benchmark NAME SOURCES)
    add_executable(Benchmark-${NAME} ${SOURCES})
    add_test(NAME Benchmark-${NAME}
             COMMAND Benchmark-${NAME})
    target_link_libraries(Benchmark-${NAME} PhysicaCore ${GPerfTools_LIBRARY})
endfunction()
