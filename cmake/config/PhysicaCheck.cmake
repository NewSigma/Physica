if(${PHYSICA_CUDA})
    set(PHYSICA_CHECK_LIBS CUDA::cudart)
endif()

macro(physica_try_compile NAME SOURCE)
    try_compile(COMPILE_RESULT_VAR
                SOURCES ${CMAKE_HOME_DIRECTORY}/cmake/config/${SOURCE}
                CMAKE_FLAGS -DINCLUDE_DIRECTORIES=${CUDAToolkit_INCLUDE_DIRS} -DCOMPILE_DEFINITIONS=-march=native
                LINK_LIBRARIES ${PHYSICA_CHECK_LIBS}
                OUTPUT_VARIABLE COMPILE_LOG)
endmacro()

macro(physica_try_run NAME SOURCE)
    try_run(RUN_RESULT_VAR COMPILE_RESULT_VAR ${CMAKE_BINARY_DIR} ${CMAKE_HOME_DIRECTORY}/cmake/config/${SOURCE}
            CMAKE_FLAGS -DINCLUDE_DIRECTORIES=${CUDAToolkit_INCLUDE_DIRS} -DCOMPILE_DEFINITIONS=-march=native
            LINK_LIBRARIES ${PHYSICA_CHECK_LIBS}
            COMPILE_OUTPUT_VARIABLE COMPILE_LOG
            RUN_OUTPUT_STDOUT_VARIABLE RUN_OUTPUT_VAR)

    if(NOT ${COMPILE_RESULT_VAR})
        message(FATAL_ERROR "Failed to compile ${NAME}: ${COMPILE_LOG}")
    endif()

    if(${RUN_RESULT_VAR} MATCHES FAILED_TO_RUN)
        message(FATAL_ERROR "Failed to run ${NAME}")
    endif()
endmacro()
