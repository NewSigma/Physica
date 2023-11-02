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

#include "LayerBase.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

namespace Physica::Core {
    template<class ScalarType> class Linear;

    namespace Internal {
        template<class T>
        class Traits<Linear<T>> {
        public:
            using ScalarType = T;
        }
    }

    template<class ScalarType>
    class Linear : public LayerBase<Linear<ScalarType>> {
        using Base = LayerBase<Linear<ScalarType>>;
        using VectorType = typename Base::VectorType;
        using MatrixType = DenseMatrix<ScalarType, MatrixOption::Row | MatrixOption::Vector>;

        MatrixType weights;
        VectorType bias;
    public:
        Linear();
        Linear(const Linear&) = default;
        Linear(Linear&&) noexcept = default;
        ~Linear() = default;
        /* Operators */
        Linear& operator=(Linear obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] inline VectorType forward(const VectorType& x) const;
        void swap(Linear& obj) noexcept;
    };

    template<class ScalarType>
    inline typename Linear<ScalarType>::VectorType Linear<ScalarType>::forward(const VectorType& x) const {
        return weights * x + bias;
    }

    template<class ScalarType>
    void Linear<ScalarType>::swap(Linear& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        weights.swap(obj.weights);
        bias.swap(obj.bias);
    }
}
