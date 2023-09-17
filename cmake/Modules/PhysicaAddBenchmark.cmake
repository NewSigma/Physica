set(PHYSICA_BENCHMARK_LIBS PhysicaCore ${GPerfTools_LIBRARY})

function(physica_add_benchmark NAME SOURCES)
    add_executable(Benchmark-${NAME} ${SOURCES} ${ARGN})
    add_test(NAME Benchmark-${NAME}
             COMMAND Benchmark-${NAME})
    target_link_libraries(Benchmark-${NAME} ${PHYSICA_BENCHMARK_LIBS})
endfunction()
