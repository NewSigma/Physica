target_sources(Benchmark PRIVATE GEMV.cpp)

if(${PHYSICA_MKL})
    target_sources(Benchmark PRIVATE GEMV_MKL.cpp)
endif()
