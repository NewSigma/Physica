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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/UnitVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RValueMatrix.h"

namespace Physica::Core {
    template<Matrix T>
    class MatrixExp : public RValueMatrix<MatrixExp<T>> {
        using This = MatrixExp<T>;
        using Base = RValueMatrix<This>;
    public:
        using typename Base::ScalarType;
    private:
        const T& m;
    public:
        MatrixExp(const T& m_);
        MatrixExp(const This&) = delete;
        MatrixExp(This&&) noexcept = delete;
        ~MatrixExp() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<Matrix M>
        void assignTo(LValueMatrix<M>& target) const;
        [[nodiscard]] ScalarType calc(size_t, size_t) const { noImpl("calc() is low performance and should be avoided"); }
        /* Getters */
        [[nodiscard]] const T& getMatrix() const noexcept { return m; }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return m.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return m.getCol(); }
    };

    template<Matrix T>
    MatrixExp<T>::MatrixExp(const T& m_) : m(m_) {
        assert(m.getRow() == m.getCol());
    }

    template<Matrix T>
    template<Matrix M>
    void MatrixExp<T>::assignTo(LValueMatrix<M>& target) const {
        for (size_t i = 0; i < getCol(); ++i)
            target.col(i) = (*this) * UnitVector<ScalarType>(i, getRow());
    }

    template<Matrix T>
    [[nodiscard]] inline auto exp(const T& m) noexcept {
        return MatrixExp<T>(m);
    }
}

namespace Physica {
    template<Physica::Core::Matrix T>
    class Traits<MatrixExp<T>> : public Traits<T> {
    public:
        constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::AnyStorage;
    };
}

#include "MatrixExpVecProduct.h"
