add_subdirectory(IO)
add_subdirectory(Math)
## Scalar
physica_add_test(Complex Scalar/Complex.cpp)
physica_add_test(Diff Scalar/Diff.cpp)
physica_add_test(Integer Scalar/Integer.cpp)
physica_add_test(MathConst Scalar/MathConst.cpp)
physica_add_test(Real Scalar/Real.cpp)
physica_add_test(SIMD Scalar/SIMD.cpp)
## ML
physica_add_test(AP ML/AP.cpp)
## Parallel
physica_add_test(ThreadPool Parallel/ThreadPool.cpp)

add_subdirectory(Physics)
