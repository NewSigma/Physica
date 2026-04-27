target_sources(Benchmark PRIVATE GEMV.cpp)

if(${PHYSICA_MKL})
    add_subdirectory(Dispatch)

    target_sources(Benchmark PRIVATE GEMV_MKL.cpp)
endif()
