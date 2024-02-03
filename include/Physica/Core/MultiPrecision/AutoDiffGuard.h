/*
 * Copyright 2023-2024 WeiBo He.
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
        constexpr static bool isDeviceSide = Utils::is_device_obj<ScalarType>::value;
        using This = AutoDiffGuard<ScalarType>;
        using PlainScalar = typename ScalarType::PlainScalar;
        using TracerType = typename std::conditional<isDeviceSide, device_obj<DiffTracer<PlainScalar>>, DiffTracer<PlainScalar>>::type;

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
        if constexpr (isDeviceSide) {
            auto& segment = TracerType::getInstance().pushSegment(ExpressionType::Set);
            node = ScalarType(segment.getValues().data(), segment.getGrads().data());
        }
        else {
            const ScalarType anyNewNode = ScalarType(0);
            node = anyNewNode;
        }
    }

    template<class ScalarType>
    AutoDiffGuard<ScalarType>::~AutoDiffGuard() {
        TracerType::getInstance().forget(node);
    }
}
