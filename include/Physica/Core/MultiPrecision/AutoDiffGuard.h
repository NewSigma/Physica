/*
 * Copyright 2023-2024 Weibo He.
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

#include "Diff.h"

namespace Physica::Core {
    template<class ScalarType>
    class AutoDiffGuard {
        static_assert(ScalarType::isDifferentiable, "[Error]: ScalarType must be differentiable");
        constexpr static bool isDeviceSide = is_device_obj<ScalarType>::value;
        using This = AutoDiffGuard<ScalarType>;
    public:
        using ValueType = typename ScalarType::ValueType;
        using TracerType = typename ScalarType::TracerType;
    private:
        ScalarType node;
    public:
        AutoDiffGuard();
        AutoDiffGuard(const AutoDiffGuard&) = delete;
        AutoDiffGuard(AutoDiffGuard&& other) noexcept;
        ~AutoDiffGuard();
        /* Operators */
        AutoDiffGuard& operator=(AutoDiffGuard obj) noexcept { swap(obj); return *this; }
        This& operator=(std::nullptr_t) noexcept { node = nullptr; return *this; }
        [[nodiscard]] operator const ScalarType&() const noexcept { return node; }
        /* Operations */
        void swap(AutoDiffGuard& __restrict obj) noexcept { node.swap(obj.node); }
        /* Getters */
        [[nodiscard]] bool isValid() const noexcept { return node != nullptr; }
    };

    template<class ScalarType>
    AutoDiffGuard<ScalarType>::AutoDiffGuard() {
        if constexpr (isDeviceSide) {
            const auto& segment = TracerType::getInstance().pushSegment(1, ExprType::Set);
            node = segment[0];
        }
        else {
            TracerType::getInstance().pushSegment(1);
            const ScalarType anyNewNode = ScalarType(0);
            node = anyNewNode;
        }
    }

    template<class ScalarType>
    AutoDiffGuard<ScalarType>::AutoDiffGuard(AutoDiffGuard&& other) noexcept : node(other.node) {
        other = nullptr;
    }

    template<class ScalarType>
    AutoDiffGuard<ScalarType>::~AutoDiffGuard() {
        if (isValid()) {
            TracerType::getInstance().forget_to(node);
            node = nullptr;
        }
    }
}
