add_subdirectory(BLAS1)
add_subdirectory(BLAS2)
add_subdirectory(BLAS3)
add_subdirectory(Eigen)
add_subdirectory(Matrix)
add_subdirectory(MatrixDecomp)
add_subdirectory(VectorExpr)

target_sources(Benchmark PRIVATE LnSumExp.cpp)
