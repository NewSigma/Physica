# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2025-2026 Weibo He
# SPDX-License-Identifier: MPL-2.0
if (${PHYSICA_HDF5})
    add_subdirectory(HDF5)
endif()

target_sources(PhysicaCore PRIVATE
    Gnuplot.cpp
    Mnist.cpp
)
