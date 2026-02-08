find_package(benchmark REQUIRED)

add_executable(Benchmark)

add_subdirectory(ManyBody)
add_subdirectory(Math)
# RPMD
target_sources(Benchmark PRIVATE RPMD/HardCore.cpp)
target_sources(Benchmark PRIVATE RPMD/ParaH.cpp)
target_sources(Benchmark PRIVATE RPMD/Q_TIP4P.cpp)
if(${PHYSICA_CUDA})
    target_sources(Benchmark PRIVATE RPMD/ParaH.cu)
endif()
# Parallel
target_sources(Benchmark PRIVATE Parallel/ThreadPool.cpp)
# Anonymous
target_sources(Benchmark PRIVATE Benchmark.cpp)
target_sources(Benchmark PRIVATE CRCoro.cpp)
target_sources(Benchmark PRIVATE Mnist.cpp)

get_target_property(benchmark_INCLUDE_DIR benchmark::benchmark INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(Benchmark PRIVATE .)
target_include_directories(Benchmark PRIVATE ${benchmark_INCLUDE_DIR}) # Explicit include to help clangd

target_link_libraries(Benchmark PhysicaCore benchmark::benchmark)
