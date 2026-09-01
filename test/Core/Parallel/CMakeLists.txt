physica_add_test(ThreadPool ThreadPool.cpp)
physica_add_test(Task Task.cpp)

if(${PHYSICA_CUDA})
    physica_add_test(ThreadBlock ThreadBlock.cu)
endif()

if(${PHYSICA_MPI})
    physica_add_mpi_test(MPI MPI.cpp 2)
endif()
