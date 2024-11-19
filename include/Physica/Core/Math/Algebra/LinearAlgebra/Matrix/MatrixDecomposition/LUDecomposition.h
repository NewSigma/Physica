/*
 * Copyright 2021-2024 Weibo He.
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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h>

namespace Physica::Core {
    template<Scalar T>
    class LUDecomposition {
        using This = LUDecomposition;
        using WorkingMatrix = DenseMatrix<T>;
    private:
        WorkingMatrix working;
    public:
        LUDecomposition(size_t order);
        template<Matrix M>
        LUDecomposition(const M& source);
        LUDecomposition(const This&) = default;
        LUDecomposition(This&&) noexcept = default;
        ~LUDecomposition() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<Matrix M>
        void compute(const M& source);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getOrder() const noexcept { return working.getRow(); }
        [[nodiscard]] size_t getRow() const noexcept { return getOrder(); }
        [[nodiscard]] size_t getCol() const noexcept { return getOrder(); }
        [[nodiscard]] const WorkingMatrix& getMatrixLU() const noexcept { return working; }
    private:
        void decompositionColumn(size_t col);
    };

    template<Scalar T>
    LUDecomposition<T>::LUDecomposition(size_t order) : working(order, order) {}

    template<Scalar T>
    template<Matrix M>
    LUDecomposition<T>::LUDecomposition(const M& source) {
        compute(source);
    }

    template<Scalar T>
    template<Matrix M>
    void LUDecomposition<T>::compute(const M& source) {
        assert(source.getRow() == source.getCol());
        working = source;
        const size_t order = getOrder();
        for (size_t i = 0; i < order; ++i)
            decompositionColumn(i);
    }

    template<Scalar T>
    void LUDecomposition<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        working.swap(obj);
    }
    /**
     * Apply LU Decomposition on a column of Matrix \from, save the result to Matrix \to.
     *
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:32
     */
    template<Scalar T>
    void LUDecomposition<T>::decompositionColumn(size_t col) {
        const auto startAlphaIndex = col + 1;
        for (size_t j = 1; j < startAlphaIndex; ++j) {
            T temp(working(j, col));
            for (size_t k = 0; k < j; ++k)
                temp -= working(j, k) * working(k, col);
            working(j, col) = std::move(temp);
        }

        const auto r = working.getRow();
        for (size_t j = startAlphaIndex; j < r; ++j) {
            T temp(working(j, col));
            for (size_t k = 0; k < col; ++k)
                temp -= working(j, k) * working(k, col);
            working(j, col) = temp / working(col, col);
        }
    }
}
