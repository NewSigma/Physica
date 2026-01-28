physica_add_test(EigenSolver Eigen/EigenSolver.cpp)
physica_add_test(JacobiDavidson Eigen/JacobiDavidson.cpp)
physica_add_test(SymmEigenSolver Eigen/SymmEigenSolver.cpp)

physica_add_test(LinearSystem LinearSystem/LinearSystem.cpp)

physica_add_test(TensorBlock Tensor/TensorBlock.cpp)

add_subdirectory(Matrix)
add_subdirectory(MatrixDecomp)

physica_add_test(Covariance Vector/Covariance.cpp)
physica_add_test(DenseVector Vector/DenseVector.cpp)
physica_add_test(DiffVector Vector/DiffVector.cpp)
physica_add_test(Givens Vector/Givens.cpp)
physica_add_test(HouseHolder Vector/HouseHolder.cpp)
if(${PHYSICA_CUDA})
    physica_add_test(DenseVector_cuda Vector/DenseVector.cu)
    physica_add_test(DiffVector_cuda Vector/DiffVector.cu)
    physica_add_test(RValueVector_cuda Vector/RValueVector.cu)
endif()
