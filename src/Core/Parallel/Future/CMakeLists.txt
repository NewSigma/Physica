target_sources(PhysicaCore PRIVATE ProcessFuture.cpp)

if(${PHYSICA_CUDA})
    target_sources(PhysicaCore PRIVATE StreamFuture.cu)
endif()
