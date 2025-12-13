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

#include "../RValueMatrix.cuh"

namespace Physica {
    template<Matrix M1, Matrix M2>
    class device_obj<Kronecker<M1, M2>> : public device_obj<RValueMatrix<Kronecker<M1, M2>>> {
        using host_obj = Kronecker<M1, M2>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
    protected:
        using typename Base::T;
    private:
        PlainStruct<const M1> m1;
        PlainStruct<const M2> m2;
    public:
        __host__ __device__ device_obj(const M1& m1, const M2& m2);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto&& target) const;
        void assign_add(Matrix auto&& target) const;

        [[nodiscard]] __device__ T calc(size_t row, size_t col) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept;
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept;
        [[nodiscard]] __host__ __device__ const auto& getLHS() const noexcept { return m1.getDerived(); }
        [[nodiscard]] __host__ __device__ const auto& getRHS() const noexcept { return m2.getDerived(); }
    };

    template<Matrix M1, Matrix M2>
    __host__ __device__ device_obj<Kronecker<M1, M2>>::device_obj(const M1& m1, const M2& m2) : m1(asStruct(m1)), m2(asStruct(m2)) {}

    template<Matrix M1, Matrix M2>
    void device_obj<Kronecker<M1, M2>>::assign(Matrix auto&& target) const {
        target.zeros();

        const auto& lhs = getLHS();
        const auto& rhs = getRHS();
        for (size_t r = 0; r < lhs.getRow(); ++r) {
            size_t offsetR = r * rhs.getRow();
            if constexpr (instanceof_tx<UnitMatrix, M1>)
                rhs.assign(target.block(offsetR, rhs.getRow(), offsetR, rhs.getCol()));
            else if constexpr (instanceof_tx<DiagMatrix, M2>)
                (rhs * lhs.calc(r, r)).assign(target.block(offsetR, rhs.getRow(), offsetR, rhs.getCol()));
            else {
                for (size_t c = 0; c < lhs.getCol(); ++c) {
                    size_t offsetC = c * rhs.getCol();
                    (rhs * lhs.calc(r, c)).assign(target.block(offsetR, rhs.getRow(), offsetC, rhs.getCol()));
                }
            }
        }
    }

    template<Matrix M1, Matrix M2>
    void device_obj<Kronecker<M1, M2>>::assign_add(Matrix auto&& target) const {
        const auto& lhs = getLHS();
        const auto& rhs = getRHS();
        for (size_t r = 0; r < lhs.getRow(); ++r) {
            size_t offsetR = r * rhs.getRow();
            if constexpr (instanceof_tx<UnitMatrix, M1>)
                rhs.assign_add(target.block(offsetR, rhs.getRow(), offsetR, rhs.getCol()));
            else if constexpr (instanceof_tx<DiagMatrix, M2>)
                (rhs * lhs.calc(r, r)).assign_add(target.block(offsetR, rhs.getRow(), offsetR, rhs.getCol()));
            else {
                for (size_t c = 0; c < lhs.getCol(); ++c) {
                    size_t offsetC = c * rhs.getCol();
                    (rhs * lhs.calc(r, c)).assign_add(target.block(offsetR, rhs.getRow(), offsetC, rhs.getCol()));
                }
            }
        }
    }

    template<Matrix M1, Matrix M2>
    __device__ auto device_obj<Kronecker<M1, M2>>::calc(size_t row, size_t col) const -> T {
        size_t row1 = row / getRHS().getRow();
        size_t row2 = row % getRHS().getRow();
        size_t col1 = col / getRHS().getCol();
        size_t col2 = col % getRHS().getCol();
        return getLHS().calc(row1, col1) * getRHS().calc(row2, col2);
    }

    template<Matrix M1, Matrix M2>
    __host__ __device__ size_t device_obj<Kronecker<M1, M2>>::getRow() const noexcept {
        return getLHS().getRow() * getRHS().getRow();
    }
    
    template<Matrix M1, Matrix M2>
    __host__ __device__ size_t device_obj<Kronecker<M1, M2>>::getCol() const noexcept {
        return getLHS().getCol() * getRHS().getCol();
    }

    template<Matrix M1, Matrix M2>
    [[nodiscard]] __host__ __device__ auto kronecker(const M1& m1, const M2& m2) noexcept requires(CUDA<M1> && CUDA<M2>) {
        return device_obj<Kronecker<M1, M2>>(m1, m2);
    }
}

namespace Physica {
    template<Matrix M1, Matrix M2>
    class Traits<device_obj<Kronecker<M1, M2>>> : public Traits<Kronecker<M1, M2>> {};
}
