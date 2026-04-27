target_sources(Benchmark PRIVATE DenseLU.cpp)
target_sources(Benchmark PRIVATE InverseLU.cpp)
target_sources(Benchmark PRIVATE DenseQR.cpp)

if(${PHYSICA_MKL})
    add_subdirectory(Dispatch)
    add_subdirectory(MKL)
endif()
