target_sources(Benchmark PRIVATE GEMV.cpp)

if(${PHYSICA_CUDA})
    target_sources(Benchmark PRIVATE GEMV.cu)
endif()
