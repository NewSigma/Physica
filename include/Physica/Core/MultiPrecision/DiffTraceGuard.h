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

#include "DifferentiableImpl/DiffTracer.h"

namespace Physica::Core {
    template<class ScalarType>
    class DiffTraceGuard {
        size_t recordIndex;
    public:
        DiffTraceGuard(const DiffTraceGuard&) = delete;
        DiffTraceGuard(DiffTraceGuard&&) noexcept = default;
        ~DiffTraceGuard();
        /* Operators */
        DiffTraceGuard& operator=(DiffTraceGuard obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(DiffTraceGuard& obj) noexcept { std::swap(recordIndex, obj.recordIndex); }
        /* Static members */
        [[nodiscard]] static DiffTraceGuard make_guard() noexcept;
    private:
        DiffTraceGuard(size_t index);
    };

    template<class ScalarType>
    DiffTraceGuard<ScalarType>::DiffTraceGuard(size_t index) : recordIndex(index) {}

    template<class ScalarType>
    DiffTraceGuard<ScalarType>::~DiffTraceGuard() {
        DiffTracer<ScalarType>::getInstance().forget(recordIndex);
    }

    template<class ScalarType>
    DiffTraceGuard<ScalarType> DiffTraceGuard<ScalarType>::make_guard() noexcept {
        return DiffTraceGuard(DiffTracer<ScalarType>::getInstance().getNumRecord());
    }
}
