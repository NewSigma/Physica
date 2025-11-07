target_sources(Benchmark PRIVATE QRDecomp.cpp)

if(${PHYSICA_MKL})
    target_sources(Benchmark PRIVATE QRDecomp_MKL.cpp)
endif()
