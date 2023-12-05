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
    template<class PlainScalar>
    class AutoDiffGuard {
        static_assert(!PlainScalar::isDifferentiable, "[Error]: Differentiable<> pack is not necessary");
        using ScalarType = Differentiable<PlainScalar, DiffMode::Reverse>;

        ScalarType node;
    public:
        AutoDiffGuard();
        AutoDiffGuard(const AutoDiffGuard&) = delete;
        AutoDiffGuard(AutoDiffGuard&&) noexcept = default;
        ~AutoDiffGuard();
        /* Operators */
        AutoDiffGuard& operator=(AutoDiffGuard obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(AutoDiffGuard& __restrict obj) noexcept { std::swap(node, obj.node); }
        /* Getters */
        [[nodiscard]] ScalarType getNode() const noexcept { return node; }
    };

    template<class PlainScalar>
    AutoDiffGuard<PlainScalar>::AutoDiffGuard() {
        const ScalarType anyNewNode = ScalarType(0);
        node = anyNewNode;
    }

    template<class PlainScalar>
    AutoDiffGuard<PlainScalar>::~AutoDiffGuard() {
        DiffTracer<PlainScalar>::getInstance().forget(node);
    }
}
