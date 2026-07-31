if(${PHYSICA_MKL}) # FIXME: Remove MKL dependency
    physica_add_test(EigenSolver Eigen/EigenSolver.cpp)
    physica_add_test(JacobiDavidson Eigen/JacobiDavidson.cpp)
endif()
physica_add_test(SymmEigenSolver Eigen/SymmEigenSolver.cpp)

physica_add_test(TensorBlock Tensor/TensorBlock.cpp)

add_subdirectory(Matrix)
add_subdirectory(MatrixDecomp)
add_subdirectory(Vector)

physica_add_test(GMRES GMRES.cpp)

if(NOT ${CMAKE_BUILD_TYPE} MATCHES "Debug")
    include(SlowTests.cmake)
endif()
