# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2024-2026 Weibo He
# SPDX-License-Identifier: GPL-3.0-or-later
physica_add_test(Ewald ForceModel/Ewald/Ewald.cpp)
physica_add_test(RandomBatchEwald ForceModel/Ewald/RandomBatchEwald.cpp)
physica_add_test(Q_TIP4P_MD ForceModel/Q_TIP4P_MD.cpp)
physica_add_test(Langevin Thermostat/Langevin.cpp)
physica_add_test(RPMD_NVE RPMD_NVE.cpp)
physica_add_test(RPMD RPMD.cpp)
if(${PHYSICA_CUDA})
    physica_add_test(RPMD_cuda RPMD.cu)
endif()
