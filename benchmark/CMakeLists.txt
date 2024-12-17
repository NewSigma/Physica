find_package(benchmark REQUIRED)

add_executable(Benchmark)

add_subdirectory(ManyBody)
# Math
target_sources(Benchmark PRIVATE Math/InnerDot.cpp)
target_sources(Benchmark PRIVATE Math/InnerDot_MKL.cpp)
target_sources(Benchmark PRIVATE Math/LnCosh.cpp)
# RPMD
target_sources(Benchmark PRIVATE RPMD/HardCore.cpp)
target_sources(Benchmark PRIVATE RPMD/ParaH.cpp)
target_sources(Benchmark PRIVATE RPMD/Q_TIP4P.cpp)
if(${PHYSICA_CUDA})
    target_sources(Benchmark PRIVATE RPMD/ParaH.cu)
    target_sources(Benchmark PRIVATE RPMD/ParaH_auto.cu)
    target_sources(Benchmark PRIVATE RPMD/Q_TIP4P2D.cu)
endif()
# Anonymous
#target_sources(Benchmark PRIVATE Mnist.cpp)

target_precompile_headers(Benchmark REUSE_FROM PhysicaCore)
get_target_property(benchmark_INCLUDE_DIR benchmark::benchmark INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(Benchmark PRIVATE ${benchmark_INCLUDE_DIR}) # Explicit include to help clangd

target_link_libraries(Benchmark PhysicaCore benchmark::benchmark_main)
