/*
 * Copyright 2023 WeiBo He.
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

#include "Differentiable.h"

namespace Physica::Core {
    template<class ScalarType>
    class AutoDiffGuard {
        static_assert(ScalarType::isDifferentiable, "[Error]: ScalarType must be differentiable");
        static_assert(!Utils::is_device_obj<ScalarType>::value, "[Error]: Include AutoDiffGuard.cuh to use the diff guard for CUDA");
        using This = AutoDiffGuard<ScalarType>;

        ScalarType node;
    public:
        AutoDiffGuard();
        AutoDiffGuard(const AutoDiffGuard&) = delete;
        AutoDiffGuard(AutoDiffGuard&&) noexcept = default;
        ~AutoDiffGuard();
        /* Operators */
        AutoDiffGuard& operator=(AutoDiffGuard obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(AutoDiffGuard& __restrict obj) noexcept { node.swap(obj.node); }
        /* Getters */
        [[nodiscard]] ScalarType getNode() const noexcept { return node; }
    };

    template<class ScalarType>
    AutoDiffGuard<ScalarType>::AutoDiffGuard() {
        const ScalarType anyNewNode = ScalarType(0);
        node = anyNewNode;
    }

    template<class ScalarType>
    AutoDiffGuard<ScalarType>::~AutoDiffGuard() {
        using PlainScalar = typename ScalarType::PlainScalar;
        DiffTracer<PlainScalar>::getInstance().forget(node);
    }
}
