if (${PHYSICA_CUDA})
    add_subdirectory(CUDA)
endif()

if (${PHYSICA_MKL})
    add_subdirectory(MKL)
endif()

if (${PHYSICA_MPI})
    target_sources(PhysicaCore PRIVATE MPIException.cpp)
endif()
