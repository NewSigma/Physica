target_sources(Benchmark PRIVATE EigenSolver.cpp)
target_sources(Benchmark PRIVATE SymmEigenSolver.cpp)

if(${PHYSICA_MKL} AND NOT ${PHYSICA_EmitLLVM})
    add_subdirectory(Dispatch)
    add_subdirectory(MKL)
endif()
