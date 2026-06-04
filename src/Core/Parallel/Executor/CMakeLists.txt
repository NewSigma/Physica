target_sources(PhysicaCore PRIVATE ProcessExecutor.cpp)
if(${PHYSICA_CUDA})
    target_sources(PhysicaCore PRIVATE CUDAExecutor.cu)
endif()
