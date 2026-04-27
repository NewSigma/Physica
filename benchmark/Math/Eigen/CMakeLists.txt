target_sources(Benchmark PRIVATE EigenSolver.cpp)
target_sources(Benchmark PRIVATE SymmEigenSolver.cpp)

if(${PHYSICA_MKL})
    add_subdirectory(Dispatch)
    add_subdirectory(MKL)
endif()
