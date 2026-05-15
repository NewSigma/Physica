physica_add_test(MatrixExp MatrixFunc/MatrixExp.cpp)
physica_add_test(MatrixPow MatrixFunc/MatrixPow.cpp)

add_subdirectory(MatrixImpl)

physica_add_test(HermiteToeplitzMatrix Toeplitz/HermiteToeplitzMatrix.cpp)

physica_add_test(DenseMatrix DenseMatrix.cpp)
physica_add_test(DenseSymmMatrix DenseSymmMatrix.cpp)
physica_add_test(DiagMatrix DiagMatrix.cpp)
physica_add_test(DiffDenseMatrix DiffDenseMatrix.cpp)
physica_add_test(MatrixChain MatrixChain.cpp)
physica_add_test(PermMatrix PermMatrix.cpp)
physica_add_test(SparseMatrix SparseMatrix.cpp)

if(${PHYSICA_CUDA})
    physica_add_test(DenseMatrix_cuda DenseMatrix.cu)
endif()
