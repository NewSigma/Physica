physica_add_test(ActionMatrix ActionMatrix.cpp)
physica_add_test(AndersonMatrix AndersonMatrix.cpp)
if(${PHYSICA_MKL}) # FIXME: Remove MKL dependency
    physica_add_test(DQMC DQMC.cpp)
    physica_add_test(FreqDQMC FreqDQMC.cpp)
endif()
physica_add_test(HubbardMatrix HubbardMatrix.cpp)
physica_add_test(J1J2Matrix J1J2Matrix.cpp)
if(${PHYSICA_MKL}) # FIXME: Remove MKL dependency
    physica_add_test(KFermiRepr KFermiRepr.cpp)
endif()
physica_add_test(SpinState SpinState.cpp)
physica_add_test(TransIsingMatrix TransIsingMatrix.cpp)

if(${PHYSICA_CUDA})
    physica_add_test(ActionMatrix_cuda ActionMatrix.cu)
    physica_add_test(FreqDQMC_cuda FreqDQMC.cu)
endif()
