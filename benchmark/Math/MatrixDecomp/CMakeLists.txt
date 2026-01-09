target_sources(Benchmark PRIVATE DenseLU/DenseLU.cpp)
target_sources(Benchmark PRIVATE DenseLU/Inverse.cpp)
target_sources(Benchmark PRIVATE DenseQR.cpp)

if(${PHYSICA_MKL})
    target_sources(Benchmark PRIVATE DenseLU/DenseLU_MKL.cpp)
    target_sources(Benchmark PRIVATE DenseLU/Inverse_MKL.cpp)
    target_sources(Benchmark PRIVATE DenseQR_MKL.cpp)
endif()
