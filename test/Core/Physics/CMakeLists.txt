# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2024-2026 Weibo He
# SPDX-License-Identifier: MPL-2.0
physica_add_test(GaussBase GaussBase.cpp)
physica_add_test(HamiltonMC HamiltonMC.cpp)
physica_add_test(IceGenerator IceGenerator.cpp)
physica_add_test(RHFSolver RHFSolver.cpp)
### DFT
physica_add_test(KSSolver ElectronicStructure/DFT/KSSolver.cpp)

add_subdirectory(ManyBody)
add_subdirectory(MD)
### Phonon
physica_add_test(FCSwapVector Phonon/FCSwapVector.cpp)
physica_add_test(PhononSolver Phonon/PhononSolver.cpp)
if(${PHYSICA_HDF5})
    physica_add_test(FrozenPhonon Phonon/FrozenPhonon.cpp)
endif()
### SolidState
physica_add_test(BrillouInterpolate SolidState/BrillouInterpolate.cpp)
physica_add_test(CrystalCell SolidState/CrystalCell.cpp)
physica_add_test(Outcar SolidState/VASP/Outcar.cpp)
physica_add_test(Poscar SolidState/VASP/Poscar.cpp)
physica_add_test(PWscfIn SolidState/QE/PWscfIn.cpp)

if(NOT ${CMAKE_BUILD_TYPE} MATCHES "Debug")
    include(SlowTests.cmake)
endif()
