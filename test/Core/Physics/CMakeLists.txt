## Physics
physica_add_test(GaussBase GaussBase.cpp)
physica_add_test(IceGenerator IceGenerator.cpp)
physica_add_test(RHFSolver RHFSolver.cpp)
### ManyBody
physica_add_test(Hubbard ManyBody/Hubbard.cpp)
physica_add_test(KSpinRepr ManyBody/KSpinRepr.cpp)
physica_add_test(SpinlessFermion ManyBody/SpinlessFermion.cpp)
### MD
physica_add_test(RSpaceEwald MD/ForceModel/Ewald/RSpaceEwald.cpp)
physica_add_test(BKSModel MD/ForceModel/BKSModel.cpp)
physica_add_test(PairModel MD/ForceModel/PairModel.cpp)
physica_add_test(SilveraGoldman MD/ForceModel/SilveraGoldman.cpp)
physica_add_test(TodaModel MD/ForceModel/TodaModel.cpp)
physica_add_test(MDCell MD/MDCell.cpp)
physica_add_test(RPMD_NVE MD/RPMD_NVE.cpp)
if(${PHYSICA_CUDA})
    physica_add_test(RSpaceEwald_cuda MD/ForceModel/Ewald/RSpaceEwald.cu)
    physica_add_test(SilveraGoldman_cuda MD/ForceModel/SilveraGoldman.cu)
    physica_add_test(HardCore_cuda MD/HardCore.cu)
endif()
### Phonon
physica_add_test(FCSwapVector Phonon/FCSwapVector.cpp)
### SolidState
physica_add_test(BrillouInterpolate SolidState/BrillouInterpolate.cpp)
physica_add_test(CrystalCell SolidState/CrystalCell.cpp)

if(NOT ${CMAKE_BUILD_TYPE} MATCHES "Debug")
    include(SlowTests.cmake)
endif()
