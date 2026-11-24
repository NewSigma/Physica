/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
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
        using ValueT = decltype([]() {
            using U = std::remove_cvref_t<Target>;
            using Tv = U::ScalarType::ValueType;
            static_assert(std::same_as<T, Tv>, "[Error]: Optimizer-Target ScalarType mismatch");
            if constexpr (Scalar<Target>)
                return T{};
            else {
                if constexpr (Vector<Target>)
                    static_assert(U::isLValueVector(), "[Error]: Cannot optimize a rvalue object");
                else if constexpr (Matrix<Target>)
                    static_assert(U::isLValueMatrix(), "[Error]: Cannot optimize a rvalue object");
                else
                    static_assert(U::isLValueTensor(), "[Error]: Cannot optimize a rvalue object");
                return typename U::template rebind_scalar<T>{};
            }
        }());
    };
}
