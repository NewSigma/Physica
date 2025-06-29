/*
 * Copyright 2023-2025 Weibo He.
 *
 * This file is part of Physica.
 *
 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include "FreeModel.h"

namespace Physica {
    template<Scalar T,
             unsigned int Dim,
             size_t NumReplica,
             RPMDIntegrator Integrator = RPMDIntegrator::Exact>
    class OpenModel : public FreeModel<T, Dim, NumReplica, Integrator> {
        using Base = FreeModel<T, Dim, NumReplica, Integrator>;
    public:
        using Base::Base;
    };

    template<Scalar T, unsigned int D, size_t N, RPMDIntegrator I>
    class Traits<OpenModel<T, D, N, I>> : public Traits<FreeModel<T, D, N, I>> {
    public:
        constexpr static bool IsPeriodBoundary = false;
    };
}
