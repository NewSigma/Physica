if (${PHYSICA_CUDA})
    target_sources(PhysicaCore PRIVATE CUDA/cuRAND.cu)
endif()

if (${PHYSICA_MKL})
    target_sources(PhysicaCore PRIVATE MKL/VSL.cpp)
endif()

if (${PHYSICA_MPI})
    target_sources(PhysicaCore PRIVATE MPIException.cpp)
endif()
