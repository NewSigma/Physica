target_sources(Benchmark PRIVATE GEMM.cpp)

if(${PHYSICA_MKL} AND NOT ${PHYSICA_EmitLLVM})
    add_subdirectory(Dispatch)

    target_sources(Benchmark PRIVATE GEMM_MKL.cpp)
endif()
