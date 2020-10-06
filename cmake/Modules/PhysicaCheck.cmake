macro(physica_check NAME)
    message(STATUS "Performing ${NAME}")
    try_run(RUN_RESULT_VAR COMPILE_RESULT_VAR ${CMAKE_BINARY_DIR} ${CMAKE_HOME_DIRECTORY}/config/Check/${NAME}.cpp)
    if(NOT ${COMPILE_RESULT_VAR})
        message(FATAL_ERROR Failed to compile ${NAME}.cpp)
    endif()

    if(${RUN_RESULT_VAR} MATCHES FAILED_TO_RUN)
        message(FATAL_ERROR Failed to run ${NAME})
    endif()
endmacro()