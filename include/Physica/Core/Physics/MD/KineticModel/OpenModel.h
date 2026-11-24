/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
        /* Static members */
        [[nodiscard]] consteval static bool isPeriodBoundary() noexcept { return false; }
    };

    template<Scalar T, unsigned int D, size_t N, RPMDIntegrator I>
    class Traits<OpenModel<T, D, N, I>> : public Traits<FreeModel<T, D, N, I>> {};
}
