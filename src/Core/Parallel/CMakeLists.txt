add_subdirectory(Executor)
add_subdirectory(Future)

target_sources(PhysicaCore PRIVATE
    MPIContext.cpp
    PBSWrapper.cpp
    SubProcess.cpp
    ThreadPool.cpp
)

if(${PHYSICA_CUDA})
    target_sources(PhysicaCore PRIVATE
        CUDAContext.cu
        CUDAEvent.cu
        CUDAStream.cu
    )
endif()
