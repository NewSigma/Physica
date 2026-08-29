add_subdirectory(Hamilton)

physica_add_test(ActionMatrix ActionMatrix.cpp)
physica_add_test(BayesDoping BayesDoping.cpp)
physica_add_test(DQMC DQMC.cpp)
physica_add_test(FreqDQMC FreqDQMC.cpp)
physica_add_test(KFermiRepr KFermiRepr.cpp)
physica_add_test(MatrixSampler MatrixSampler.cpp)
physica_add_test(SpinState SpinState.cpp)

if(${PHYSICA_CUDA})
    physica_add_test(ActionMatrix_cuda ActionMatrix.cu)
    physica_add_test(FreqDQMC_cuda FreqDQMC.cu)
endif()
