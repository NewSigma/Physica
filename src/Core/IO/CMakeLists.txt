# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2025-2026 Weibo He
# SPDX-License-Identifier: GPL-3.0-or-later
if (${PHYSICA_HDF5})
    add_subdirectory(HDF5)
endif()

target_sources(PhysicaCore PRIVATE
    Gnuplot.cpp
    Mnist.cpp
)
