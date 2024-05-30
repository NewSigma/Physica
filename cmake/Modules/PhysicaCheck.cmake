if(${PHYSICA_CUDA})
    set(PHYSICA_CHECK_LIBS CUDA::cudart)
endif()

macro(physica_check NAME SOURCE)
    try_run(RUN_RESULT_VAR COMPILE_RESULT_VAR ${CMAKE_BINARY_DIR} ${CMAKE_HOME_DIRECTORY}/cmake/config/${SOURCE} 
            CMAKE_FLAGS -DINCLUDE_DIRECTORIES=${CUDAToolkit_INCLUDE_DIRS}
            LINK_LIBRARIES ${PHYSICA_CHECK_LIBS})

    if(NOT ${COMPILE_RESULT_VAR})
        message(FATAL_ERROR "Failed to compile ${NAME}")
    endif()

    if(${RUN_RESULT_VAR} MATCHES FAILED_TO_RUN)
        message(FATAL_ERROR "Failed to run ${NAME}")
    endif()
endmacro()