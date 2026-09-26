physica_add_test(EigenSolver Eigen/EigenSolver.cpp)
physica_add_test(JacobiDavidson Eigen/JacobiDavidson.cpp)
physica_add_test(SymmEigenSolver Eigen/SymmEigenSolver.cpp)

physica_add_test(DenseTensor Tensor/DenseTensor.cpp)
physica_add_test(LValueTensor Tensor/LValueTensor.cpp)
physica_add_test(CompactTensor Tensor/CompactTensor.cpp)
physica_add_test(RValueTensor Tensor/RValueTensor.cpp)

add_subdirectory(Matrix)
add_subdirectory(MatrixDecomp)
add_subdirectory(Vector)

if(${PHYSICA_CUDA})
    physica_add_test(DenseTensor_cuda Tensor/DenseTensor.cu)
endif()

physica_add_test(GMRES GMRES.cpp)

if(NOT ${CMAKE_BUILD_TYPE} MATCHES "Debug")
    include(SlowTests.cmake)
endif()
