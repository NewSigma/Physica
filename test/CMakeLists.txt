include(CTest)
add_custom_target(Tests COMMAND ctest --parallel WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

function(physica_add_test NAME SOURCES)
    add_executable(Test-${NAME})
    target_sources(Test-${NAME} PRIVATE ${SOURCES} ${ARGN})
    target_include_directories(Test-${NAME} PRIVATE ${CMAKE_SOURCE_DIR}/test)
    target_link_libraries(Test-${NAME} PhysicaCore ${GPerfTools_LIBRARY})

    add_test(NAME Test-${NAME} COMMAND Test-${NAME})

    add_dependencies(Tests Test-${NAME})
endfunction()

message(STATUS "Configuring tests.")

add_subdirectory(Core)

physica_add_test(Logger Logger/LoggerTest.cpp)
target_link_libraries(Test-Logger PhysicaLogger)
