/*
 * Copyright 2024 Weibo He.
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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RValueMatrix.h>

namespace Physica::Core {
    template<class MatrixType>
    class MatrixExp : public RValueMatrix<MatrixExp<MatrixType>> {
        using This = MatrixExp<MatrixType>;
        using Base = RValueMatrix<This>;
    public:
        using typename Base::ScalarType;
    private:
        const MatrixType& m;
    public:
        MatrixExp(const RValueMatrix<MatrixType>& m_);
        MatrixExp(const This&) = delete;
        MatrixExp(This&&) noexcept = delete;
        ~MatrixExp() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t, size_t) const { noImpl(); }
        /* Getters */
        [[nodiscard]] const MatrixType& getMatrix() const noexcept { return m; }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return m.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept { return m.getColumn(); }
    };

    template<class MatrixType>
    MatrixExp<MatrixType>::MatrixExp(const RValueMatrix<MatrixType>& m_) : m(m_.getDerived()) {
        assert(m.getRow() == m.getColumn());
    }

    template<class MatrixType>
    [[nodiscard]] inline MatrixExp<MatrixType> exp(const RValueMatrix<MatrixType>& m) noexcept {
        return MatrixExp<MatrixType>(m);
    }
}

namespace Physica {
    using namespace Core;

    template<class MatrixType>
    class Traits<MatrixExp<MatrixType>> : public Traits<MatrixType> {
    public:
        constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::AnyStorage;
    };
}

#include "MatrixExpVecProduct.h"
