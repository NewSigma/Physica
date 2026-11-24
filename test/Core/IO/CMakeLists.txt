# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2024-2026 Weibo He
# SPDX-License-Identifier: MPL-2.0
physica_add_test(Gnuplot Gnuplot.cpp)
if(${PHYSICA_HDF5})
    physica_add_test(HDF5 HDF5.cpp)
endif()
