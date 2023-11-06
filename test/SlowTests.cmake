# Core
## Math
### Calculus
physica_add_test(FokkerPlanck Core/Math/Calculus/PDE/FEM/FokkerPlanck.cpp)
### Optimization
physica_add_test(Adam Core/Math/Optimization/AdamTest.cpp)
## Physics
physica_add_test(DimEstimator Core/Physics/DimEstimator.cpp)
physica_add_test(KSSolver Core/Physics/KSSolver.cpp)
### MD
physica_add_test(Ewald Core/Physics/MD/ForceModel/Ewald.cpp)
physica_add_test(Q_TIP4P Core/Physics/MD/ForceModel/Q_TIP4P.cpp)
physica_add_test(Langevin Core/Physics/MD/Thermostat/Langevin.cpp)
physica_add_test(RPMD Core/Physics/MD/RPMD.cpp)

if(${PHYSICA_CUDA})
    physica_add_test(SilveraGoldman_cuda Core/Physics/MD/ForceModel/SilveraGoldman.cu)
    physica_add_test(RPMD_cuda Core/Physics/MD/RPMD.cu)
endif()
