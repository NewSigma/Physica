target_sources(Benchmark PRIVATE DenseLU.cpp)
target_sources(Benchmark PRIVATE DenseQR.cpp)

if(${PHYSICA_MKL})
    target_sources(Benchmark PRIVATE DenseLU_MKL.cpp)
    target_sources(Benchmark PRIVATE DenseQR_MKL.cpp)
endif()
