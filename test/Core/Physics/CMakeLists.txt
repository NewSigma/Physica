## Physics
physica_add_test(GaussBase GaussBase.cpp)
physica_add_test(IceGenerator IceGenerator.cpp)
physica_add_test(RHFSolver RHFSolver.cpp)
### ManyBody
physica_add_test(HubbardMatrix ManyBody/HubbardMatrix.cpp)
physica_add_test(KSpinRepr ManyBody/KSpinRepr.cpp)
physica_add_test(SpinlessFermion ManyBody/SpinlessFermion.cpp)
### MD
add_subdirectory(MD)
### Phonon
physica_add_test(FCSwapVector Phonon/FCSwapVector.cpp)
### SolidState
physica_add_test(BrillouInterpolate SolidState/BrillouInterpolate.cpp)
physica_add_test(CrystalCell SolidState/CrystalCell.cpp)

if(NOT ${CMAKE_BUILD_TYPE} MATCHES "Debug")
    include(SlowTests.cmake)
endif()
