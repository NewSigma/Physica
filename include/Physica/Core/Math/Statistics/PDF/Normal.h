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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

namespace Physica::Core {
    template<class ScalarType>
    class Normal {
        ScalarType mean;
        ScalarType deviation;
    public:
        Normal(ScalarType mean_ = 0, ScalarType deviation_ = 1);
        Normal(const Normal&) = default;
        Normal(Normal&&) noexcept = default;
        ~Normal() = default;
        /* Operator */
        Normal& operator=(Normal obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] ScalarType operator()(const ScalarType& x) const;
        template<class VectorType>
        [[nodiscard]] Vector<ScalarType> operator()(const RValueVector<VectorType>& x) const;
        /* Operations */
        void swap(Normal& obj) noexcept;
    };

    template<class ScalarType>
    Normal<ScalarType>::Normal(ScalarType mean_, ScalarType deviation_)
            : mean(std::move(mean_)), deviation(std::move(deviation_)) {}

    template<class ScalarType>
    ScalarType Normal<ScalarType>::operator()(const ScalarType& x) const {
        const ScalarType repDevia = reciprocal(repDevia);
        const ScalarType factor = repDevia / sqrt(ScalarType(2 * M_PI));
        return exp(square((x - mean) * repDevia) * ScalarType(-0.5)) * factor;
    }

    template<class ScalarType>
    template<class VectorType>
    Vector<ScalarType> Normal<ScalarType>::operator()(const RValueVector<VectorType>& x) const {
        const ScalarType repDevia = reciprocal(deviation);
        const ScalarType factor = repDevia / sqrt(ScalarType(2 * M_PI));
        return exp(square((x - mean) * repDevia) * ScalarType(-0.5)) * factor;
    }

    template<class ScalarType>
    void Normal<ScalarType>::swap(Normal& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        mean.swap(obj.mean);
        deviation.swap(obj.deviation);
    }
}
