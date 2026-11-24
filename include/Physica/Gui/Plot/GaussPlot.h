/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Plot.h"

namespace Physica {
    class PHYSICA_API GaussPlot : public Plot {
        using Base = Plot;
    public:
        GaussPlot(double maxX, double maxY, double deltaX, double deltaY, unsigned int numSigma);
    };
}
