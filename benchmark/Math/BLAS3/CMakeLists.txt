target_sources(Benchmark PRIVATE GEMM.cpp)

if(${PHYSICA_MKL})
    add_subdirectory(Dispatch)

    target_sources(Benchmark PRIVATE GEMM_MKL.cpp)
endif()
