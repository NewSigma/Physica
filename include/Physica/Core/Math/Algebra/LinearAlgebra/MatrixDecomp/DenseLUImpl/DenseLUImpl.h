/*
 * Copyright 2025-2026 Weibo He.
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

namespace Physica {
    template<Scalar T, bool Pivot>
    DenseLU<T, Pivot>::DenseLU(size_t size) : working(size, size), perm(size) {}

    template<Scalar T, bool Pivot>
    DenseLU<T, Pivot>::DenseLU(const Matrix auto& source) : DenseLU(source.getRow()) {
        compute(source);
    }

    template<Scalar T, bool Pivot>
    void DenseLU<T, Pivot>::compute() {
        if constexpr (HasMKL())
            compute_mkl();
        else
            compute_base();
    }

    template<Scalar T, bool Pivot>
    void DenseLU<T, Pivot>::compute_base() {
        decomp(working, 0);

        if constexpr (Pivot)
            perm = perm.inv();
    }

    template<Scalar T, bool Pivot>
    void DenseLU<T, Pivot>::compute(const Matrix auto& source) {
        source.assign(working);
        compute();
    }

    template<Scalar T, bool Pivot>
    void DenseLU<T, Pivot>::compute_base(const Matrix auto& source) {
        source.assign(working);
        compute_base();
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
    auto& DenseLU<T, Pivot>::getMatrixLU(this auto&& self) noexcept {
        return self.working;
    }

    template<Scalar T, bool Pivot>
    const auto& DenseLU<T, Pivot>::getPerm() const noexcept {
        static_assert(Pivot, "[Error]: Perm is available to PLU decomp only");
        return perm;
    }

    template<Scalar T, bool Pivot>
    void DenseLU<T, Pivot>::decomp(Matrix auto& block, size_t offset) {
        size_t order = block.getOrder();
        for (size_t i = 0; i < order; ++i) {
            auto col = block.col(i);
            for (size_t k = 0; k < i; ++k)
                col.tail(k + 1) -= block.col(k).tail(k + 1) * col[k];

            const size_t alpha = i + 1;
            if (alpha < order) {
                if constexpr (Pivot) {
                    size_t pivot = block.pivotPartial(i);
                    if (i != pivot) {
                        i += offset;
                        pivot += offset;
                        working.swap_row(i, pivot);
                        perm.swap_row(i, pivot);
                    }
                }

                col.tail(alpha) *= reciprocal(col[i]);
            }
        }
    }
}
