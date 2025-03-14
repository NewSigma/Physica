physica_add_example(DenseVector DenseVector.cpp)
physica_add_example(JacobiDavidson JacobiDavidson.cpp)

if(${PHYSICA_CUDA})
    physica_add_example(DenseVector_cuda DenseVector.cu)
endif()
