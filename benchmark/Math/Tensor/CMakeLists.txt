if(${PHYSICA_CUDA})
    target_sources(Benchmark PRIVATE StridedTensor.cu)
endif()
