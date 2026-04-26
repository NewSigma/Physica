target_sources(Benchmark PRIVATE GEMM.cpp)

if(${PHYSICA_MKL})
    target_sources(Benchmark PRIVATE GEMM_MKL.cpp)
endif()
