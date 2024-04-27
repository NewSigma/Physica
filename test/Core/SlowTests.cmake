# Core
## Math
### Optimization
physica_add_test(Adam Math/Optimization/AdamTest.cpp)
## Physics
physica_add_test(DimEstimator Physics/DimEstimator.cpp)
physica_add_test(KSSolver Physics/KSSolver.cpp)
### MD
physica_add_test(Ewald Physics/MD/ForceModel/Ewald/Ewald.cpp)
physica_add_test(RandomBatchEwald Physics/MD/ForceModel/Ewald/RandomBatchEwald.cpp)
physica_add_test(Q_TIP4P Physics/MD/ForceModel/Q_TIP4P.cpp)
physica_add_test(Q_TIP4P_MD Physics/MD/ForceModel/Q_TIP4P_MD.cpp)
physica_add_test(Langevin Physics/MD/Thermostat/Langevin.cpp)
physica_add_test(RPMD Physics/MD/RPMD.cpp)

if(${PHYSICA_CUDA})
    physica_add_test(RPMD_cuda Physics/MD/RPMD.cu)
endif()
