target_sources(Benchmark PRIVATE EigenSolver.cpp)
target_sources(Benchmark PRIVATE SymmEigenSolver.cpp)

if(${PHYSICA_MKL})
    target_sources(Benchmark PRIVATE EigenSolver_MKL.cpp)
    target_sources(Benchmark PRIVATE SymmEigenSolver_MKL.cpp)
endif()
