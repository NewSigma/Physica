/*
 * Copyright 2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/Tensor.h"

namespace Physica {
    template<Scalar T>
    class OptBase {
        static_assert(!Diffable<T>);
    protected:
        template<class Target>
        using ValueT = decltype([](const auto& target) {
            using Tv = std::remove_cvref_t<Target>::ScalarType::ValueType;
            static_assert(std::same_as<T, Tv>, "[Error]: Optimizer-Target ScalarType mismatch");
            if constexpr (Scalar<Target>)
                return T{};
            else {
                if constexpr (Vector<Target>)
                    static_assert(target.isLValueVector(), "[Error]: Cannot optimize a rvalue object");
                else if constexpr (Matrix<Target>)
                    static_assert(target.isLValueMatrix(), "[Error]: Cannot optimize a rvalue object");
                else
                    static_assert(target.isLValueTensor(), "[Error]: Cannot optimize a rvalue object");
                return typename Target::template rebind_scalar<T>{};
            }
        }());
    };
}
