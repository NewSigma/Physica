# MD
physica_add_test(Ewald MD/ForceModel/Ewald/Ewald.cpp)
physica_add_test(RandomBatchEwald MD/ForceModel/Ewald/RandomBatchEwald.cpp)
physica_add_test(Q_TIP4P MD/ForceModel/Q_TIP4P.cpp)
physica_add_test(Q_TIP4P_MD MD/ForceModel/Q_TIP4P_MD.cpp)
physica_add_test(Langevin MD/Thermostat/Langevin.cpp)
physica_add_test(RPMD MD/RPMD.cpp)

if(${PHYSICA_CUDA})
    physica_add_test(RPMD_cuda MD/RPMD.cu)
endif()
# Anonymous
physica_add_test(DimEstimator DimEstimator.cpp)
physica_add_test(KSSolver KSSolver.cpp)