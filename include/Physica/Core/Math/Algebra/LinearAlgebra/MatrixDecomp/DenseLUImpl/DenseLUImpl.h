/*
 * Copyright 2025 Weibo He.
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

#include "../DenseLU.h"
#include "MatrixL.h"

namespace Physica {
    template<Scalar T, bool Pivot>
    DenseLU<T, Pivot>::DenseLU(size_t size) : working(size, size), perm(size) {}

    template<Scalar T, bool Pivot>
    DenseLU<T, Pivot>::DenseLU(const Matrix auto& source) : DenseLU(source.getRow()) {
        compute(source);
    }

    template<Scalar T, bool Pivot>
    void DenseLU<T, Pivot>::compute(const Matrix auto& source) {
        if constexpr (HasMKL())
            compute_mkl(source);
        else
            compute_base(source);
    }

    template<Scalar T, bool Pivot>
    void DenseLU<T, Pivot>::compute_base(const Matrix auto& source) {
        source.assert_assign(working);
        source.assign(working);
        for (size_t i = 0; i < getOrder(); ++i) {
            if constexpr (Pivot) {
                size_t j = working.partialPivoting(i);
                perm.swapRows(i, j);
            }
            decomp_col(i);
        }

        if constexpr (Pivot)
            perm = perm.inv();
    }

    template<Scalar T, bool Pivot>
    T DenseLU<T, Pivot>::det() const noexcept {
        T result = working.diag().prod();
        if constexpr (Pivot)
            result *= perm.det();
        return result;
    }

    template<Scalar T, bool Pivot>
    auto DenseLU<T, Pivot>::lnAbsDet() const noexcept -> Tr {
        return ln(abs(working.diag())).sum();
    }

    template<Scalar T, bool Pivot>
    T DenseLU<T, Pivot>::sgndet() const noexcept {
        return unit(working.diag()).prod();
    }

    template<Scalar T, bool Pivot> 
    auto DenseLU<T, Pivot>::inv() const noexcept {
        return Inverse<This>(*this);
    }

    template<Scalar T, bool Pivot>
    void DenseLU<T, Pivot>::resize(size_t size) {
        working.resize(size, size);
        if constexpr (Pivot)
            perm.resize(size);
    }

    template<Scalar T, bool Pivot>
    void DenseLU<T, Pivot>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        working.swap(obj.working);
        if constexpr (Pivot)
            perm.swap(obj.perm);
    }

    template<Scalar T, bool Pivot>
    const auto& DenseLU<T, Pivot>::getPerm() const noexcept {
        static_assert(Pivot, "[Error]: Perm is available to PLU decomp only");
        return perm;
    }

    template<Scalar T, bool Pivot>
    void DenseLU<T, Pivot>::pre_compute([[maybe_unused]] const Matrix auto& source) const noexcept {
        assert(source.isSquare());
        assert(source.getRow() == getOrder());
    }
    /**
     * Reference:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005:32
     */
    template<Scalar T, bool Pivot>
    void DenseLU<T, Pivot>::decomp_col(size_t col) {
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
