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
        static_assert(!ScalarType::isDifferentiable, "[Error]: Differentiable<> pack is not necessary");

        size_t traceIndex;
    public:
        AutoDiffGuard();
        AutoDiffGuard(const AutoDiffGuard&) = delete;
        AutoDiffGuard(AutoDiffGuard&&) noexcept = default;
        ~AutoDiffGuard();
        /* Operators */
        AutoDiffGuard& operator=(AutoDiffGuard obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(AutoDiffGuard& obj) noexcept { std::swap(traceIndex, obj.traceIndex); }
        /* Getters */
        [[nodiscard]] size_t getTraceIndex() const noexcept { return traceIndex; }
    private:
        AutoDiffGuard(size_t index);
    };

    template<class ScalarType>
    AutoDiffGuard<ScalarType>::AutoDiffGuard() : traceIndex(DiffTracer<ScalarType>::getInstance().getNumRecord()) {}

    template<class ScalarType>
    AutoDiffGuard<ScalarType>::AutoDiffGuard(size_t index) : traceIndex(index) {}

    template<class ScalarType>
    AutoDiffGuard<ScalarType>::~AutoDiffGuard() {
        DiffTracer<ScalarType>::getInstance().forget(traceIndex);
    }
}
