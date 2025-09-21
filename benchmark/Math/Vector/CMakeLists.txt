target_sources(Benchmark PRIVATE Exp.cpp)
target_sources(Benchmark PRIVATE LnCosh.cpp)

if(${PHYSICA_MKL})
    target_sources(Benchmark PRIVATE Exp_MKL.cpp)
endif()
