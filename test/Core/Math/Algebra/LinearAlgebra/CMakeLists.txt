physica_add_test(EigenSolver Eigen/EigenSolver.cpp)
physica_add_test(JacobiDavidson Eigen/JacobiDavidson.cpp)
physica_add_test(SymmEigenSolver Eigen/SymmEigenSolver.cpp)

physica_add_test(TensorBlock Tensor/TensorBlock.cpp)

add_subdirectory(Matrix)
add_subdirectory(MatrixDecomp)
add_subdirectory(Vector)
