add_subdirectory(Executor)
add_subdirectory(Future)

target_sources(PhysicaCore PRIVATE
    Task/Task.cpp
    MPI.cpp
    SubProcess.cpp
    ThreadPool.cpp
)

if(${PHYSICA_CUDA})
    target_sources(PhysicaCore PRIVATE
        CUDAContext.cu
        CUDAEvent.cu
        CUDAStream.cu
        KernelConfig.cu
    )
endif()
