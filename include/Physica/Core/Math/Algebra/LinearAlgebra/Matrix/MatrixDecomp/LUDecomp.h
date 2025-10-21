/*
 * Copyright 2021-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/PermMatrix.h"

namespace Physica {
    template<Scalar T, bool Pivot>
    class LUDecomp {
        using This = LUDecomp;
        using WorkingMatrix = DenseMatrix<T, MatrixOption::Col>;
        using BiasArray = std::conditional<Pivot, PermMatrix<T>, PlainStruct<void>>::type;

        constexpr static bool isComplex = T::isComplex;
        using Tc = T::ComplexType;
        using Tm = std::conditional<isComplex, typename Tc::MKL_Complex, typename T::MachineType>::type;
    private:
        WorkingMatrix working;
        [[no_unique_address]] PermMatrix<T> perm;
    public:
        LUDecomp() = default;
        LUDecomp(size_t size);
        LUDecomp(const Matrix auto& source);
        LUDecomp(const This&) = default;
        LUDecomp(This&&) noexcept = default;
        ~LUDecomp() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void compute(const Matrix auto& source);
        void compute_mkl(const Matrix auto& source);
        void compute_base(const Matrix auto& source);

        void resize(size_t size);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getOrder() const noexcept { return working.getRow(); }
        [[nodiscard]] size_t getRow() const noexcept { return getOrder(); }
        [[nodiscard]] size_t getCol() const noexcept { return getOrder(); }
        [[nodiscard]] const auto& getMatrixLU() const noexcept { return working; }
        [[nodiscard]] auto getMatrixU() const noexcept { return working.triu(); }
        [[nodiscard]] const auto& getPerm() const noexcept { return perm; }
    private:
        void pre_compute(const Matrix auto& source) const noexcept;
        void decomp_col(size_t col);
    };

    template<Scalar T, bool Pivot>
    LUDecomp<T, Pivot>::LUDecomp(size_t size) {
        resize(size);
    }

    template<Scalar T, bool Pivot>
    LUDecomp<T, Pivot>::LUDecomp(const Matrix auto& source) : LUDecomp(source.getRow()) {
        compute(source);
    }

    template<Scalar T, bool Pivot>
    void LUDecomp<T, Pivot>::compute(const Matrix auto& source) {
        if constexpr (HasMKL())
            compute_mkl(source);
        else
            compute_base(source);
    }

    template<Scalar T, bool Pivot>
    void LUDecomp<T, Pivot>::compute_base(const Matrix auto& source) {
        pre_compute(source);

        const size_t order = source.getRow();
        if constexpr (Pivot)
            perm = PermMatrix<T>(order);

        working = source;
        for (size_t i = 0; i < order; ++i) {
            if constexpr (Pivot) {
                size_t j = working.partialPivoting(i);
                perm.swapRows(i, j);
            }
            decomp_col(i);
        }
    }

    template<Scalar T, bool Pivot>
    void LUDecomp<T, Pivot>::resize(size_t size) {
        working.resize(size, size);
        if constexpr (Pivot)
            perm.resize(size);
    }

    template<Scalar T, bool Pivot>
    void LUDecomp<T, Pivot>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        working.swap(obj);
        if constexpr (Pivot)
            perm.swap(perm);
    }

    template<Scalar T, bool Pivot>
    void LUDecomp<T, Pivot>::pre_compute(const Matrix auto& source) const noexcept {
        assert(source.isSquare());
        assert(source.getRow() == getOrder());
    }
    /**
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:32
     */
    template<Scalar T, bool Pivot>
    void LUDecomp<T, Pivot>::decomp_col(size_t col) {
        const size_t alpha = col + 1;
        for (size_t j = 1; j < alpha; ++j)
            working(j, col) -= working.row(j).head(j) * working.col(col).head(j);

        const T factor = reciprocal(working(col, col));
        if (col == 0)
            working.col(0).tail(alpha) *= factor;
        else {
            const size_t r = working.getRow();
            for (size_t j = alpha; j < r; ++j)
                working(j, col) = (working(j, col) - working.row(j).head(col) * working.col(col).head(col)) * factor;
        }
    }
}

#ifdef PHYSICA_MKL
    #include "LUDecomp_MKL.h"
#endif
