# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2026 Weibo He
# SPDX-License-Identifier: MPL-2.0
add_subdirectory(PlotImpl)

target_sources(PhysicaGui PRIVATE
    DensityPlot.cpp
    GaussPlot.cpp
    MultiPlot.cpp
    Plot3D.cpp
    Plot.cpp
    PolarPlot.cpp
)
