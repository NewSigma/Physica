# Core
## Math
### Optimization
physica_add_test(Adam Core/Math/Optimization/AdamTest.cpp)
## Physics
physica_add_test(DimEstimator Core/Physics/DimEstimator.cpp)
physica_add_test(KSSolver Core/Physics/KSSolver.cpp)
### MD
physica_add_test(Ewald Core/Physics/MD/ForceModel/Ewald/Ewald.cpp)
physica_add_test(RandomBatchEwald Core/Physics/MD/ForceModel/Ewald/RandomBatchEwald.cpp)
physica_add_test(Q_TIP4P Core/Physics/MD/ForceModel/Q_TIP4P.cpp)
physica_add_test(Q_TIP4P_MD Core/Physics/MD/ForceModel/Q_TIP4P_MD.cpp)
physica_add_test(Langevin Core/Physics/MD/Thermostat/Langevin.cpp)
physica_add_test(RPMD Core/Physics/MD/RPMD.cpp)

if(${PHYSICA_CUDA})
    physica_add_test(RPMD_cuda Core/Physics/MD/RPMD.cu)
endif()
