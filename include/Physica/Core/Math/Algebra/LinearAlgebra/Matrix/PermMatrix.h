/*
 * Copyright 2022-2026 Weibo He.
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

#include <unordered_set>
#include "MatrixImpl/RValueMatrix.h"

namespace Physica {
    template<Scalar T>
    class PermMatrix : public RValueMatrix<PermMatrix<T>> {
        using This = PermMatrix<T>;
        using Base = RValueMatrix<This>;
        static_assert(!T::isComplex(), "[Error]: Permutation matrix is real");

        using typename Base::Trv;
    private:
        Array<size_t> indices;
    public:
        PermMatrix() = default;
        PermMatrix(size_t order);
        PermMatrix(Array<size_t> indices_);
        PermMatrix(const This&) = default;
        PermMatrix(This&&) noexcept = default;
        ~PermMatrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;

        [[nodiscard]] T det() const;
        [[nodiscard]] constexpr static Trv lnAbsDet() noexcept { return Trv(0); }

        [[nodiscard]] PermMatrix inv() const noexcept;
        [[nodiscard]] PermMatrix transpose() const noexcept;
        [[nodiscard]] Array<MKL_INT64> toMKL() const;

        void resize(size_t order);
        void swap_row(size_t row1, size_t row2);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] auto&& getIndices(this auto&&) noexcept;
        [[nodiscard]] size_t getRow() const noexcept { return indices.getLength(); }
        [[nodiscard]] size_t getCol() const noexcept { return indices.getLength(); }
        /* Static members */
        [[nodiscard]] static This fromMKL(Array<MKL_INT64> ipiv);
    };

    template<Scalar T>
    PermMatrix<T>::PermMatrix(size_t order) : indices(order) {
        for (size_t i = 0; i < order; ++i)
            indices[i] = i;
    }

    template<Scalar T>
    PermMatrix<T>::PermMatrix(Array<size_t> indices_) : indices(std::move(indices_)) {
        std::unordered_set<size_t> buffer{};
        for (size_t index : indices) {
            if (index >= indices.getLength())
                throw std::invalid_argument("[Error]: Invalid index");
            buffer.insert(index);
        }
        if (buffer.size() != indices.getLength())
            throw std::invalid_argument("[Error]: Duplicate index is not allowed");
    }

    template<Scalar T>
    T PermMatrix<T>::calc(size_t row, size_t col) const {
        return indices[row] == col ? T(1) : T(0);
    }

    template<Scalar T>
    T PermMatrix<T>::det() const {
        int count = 0;
        for (size_t i = 0; i < getRow(); ++i) {
            for (size_t j = i + 1; j < getRow(); ++j) {
                count += indices[j] < indices[i];
            }
        }
        return T((count % 2 == 0) ? 1.0 : -1.0);
    }

    template<Scalar T>
    auto PermMatrix<T>::inv() const noexcept -> This {
        Array<size_t> result(indices.getLength());
        for (size_t i = 0; i < result.getLength(); ++i)
            result[indices[i]] = i;
        return This(std::move(result));
    }

    template<Scalar T>
    auto PermMatrix<T>::transpose() const noexcept -> This {
        return inv();
    }

    template<Scalar T>
    Array<MKL_INT64> PermMatrix<T>::toMKL() const {
        size_t length = indices.getLength();
        return Array<MKL_INT64>::generate(length, [perm = *this, length](MKL_INT64 i) mutable {
            for (auto j = i; j < length; ++j) {
                if (perm.getIndices()[j] == i) {
                    perm.swap_row(i, j);
                    return j + 1;
                }
            }
            unreachable();
        });
    }

    template<Scalar T>
    void PermMatrix<T>::resize(size_t order) {
        *this = PermMatrix<T>(order);
    }

    template<Scalar T>
    void PermMatrix<T>::swap_row(size_t row1, size_t row2) {
        std::swap(indices[row1], indices[row2]);
    }

    template<Scalar T>
    void PermMatrix<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        indices.swap(obj.indices);
    }

    template<Scalar T>
    auto&& PermMatrix<T>::getIndices(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).indices;
    }

    template<Scalar T>
    auto PermMatrix<T>::fromMKL(Array<MKL_INT64> ipiv) -> This {
        auto perm = This(ipiv.getLength());
        for (size_t i = 0; i < ipiv.getLength(); ++i)
            perm.swap_row(i, ipiv[i] - 1);
        return perm.inv();
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<PermMatrix<T>> {
    public:
        using ScalarType = T;
        constexpr static int Major = MatrixMajor::BothMajor;
        constexpr static size_t RowAtCompile = Dynamic;
        constexpr static size_t ColAtCompile = Dynamic;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}
