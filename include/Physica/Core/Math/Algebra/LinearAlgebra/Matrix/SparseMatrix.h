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

#include "MatrixImpl/RValueMatrix.h"

namespace Physica {
    /**
     * CSR & CSC format
     */
    template<Scalar T, int Major = MatrixMajor::Row>
    class SparseMatrix : public RValueMatrix<SparseMatrix<T, Major>> {
        using This = SparseMatrix<T, Major>;
        using Base = RValueMatrix<This>;
    private:
        Array<T> elements;
        Array<size_t> minorIndexes;
        Array<size_t> majorStarts;
        size_t maxMinor = 0;
    public:
        SparseMatrix() = default;
        SparseMatrix(size_t row, size_t col);
        explicit SparseMatrix(const Matrix auto& m);
        SparseMatrix(const This&) = default;
        SparseMatrix(This&&) noexcept = default;
        ~SparseMatrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        This& operator=(const Matrix auto& m);
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;
        void insert(T x, size_t row, size_t col);

        using Base::resize;
        void resize(size_t row, size_t col);
        [[nodiscard]] device_obj<This> toDevice() const;
        [[nodiscard]] device_obj<This> toDeviceAsync() const;
        void toDevice(device_obj<This>& obj) const;
        void toDeviceAsync(device_obj<This>& obj) const;

        void zeros();
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept;
        [[nodiscard]] size_t getCol() const noexcept;
        [[nodiscard]] size_t getMaxMajor() const noexcept;
        [[nodiscard]] size_t getMaxMinor() const noexcept;
        [[nodiscard]] size_t getNumNonZero() const noexcept { return elements.getLength(); }
        [[nodiscard]] const auto& getElements() const { return elements; }
        [[nodiscard]] const auto& getMinorIndexes() const { return minorIndexes; }
        [[nodiscard]] const auto& getMajorStarts() const { return majorStarts; }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return Major; }
        /* Friends */
        friend class device_obj<This>;
    };

    template<Scalar T, int Major>
    SparseMatrix<T, Major>::SparseMatrix(size_t row, size_t col)
            : majorStarts(MatrixMajor::selectMajor<This>(row, col) + 1, 0)
            , maxMinor(MatrixMajor::selectMinor<This>(row, col)) {
        size_t size = std::max(row, col);
        elements.reserve(size);
        minorIndexes.reserve(size);
    }

    template<Scalar T, int Major>
    SparseMatrix<T, Major>::SparseMatrix(const Matrix auto& m) : SparseMatrix(m.getRow(), m.getCol()) {
        operator=(m);
    }

    template<Scalar T, int Major>
    auto SparseMatrix<T, Major>::operator=(const Matrix auto& m) -> This& {
        resize(m);
        for (size_t r = 0; r < getRow(); ++r) {
            for (size_t c = 0; c < getCol(); ++c) {
                T x = m.calc(r, c);
                if (!x.isZero())
                    insert(x, r, c);
            }
        }
        return *this;
    }

    template<Scalar T, int Major>
    T SparseMatrix<T, Major>::calc(size_t row, size_t col) const {
        assert(row < getRow());
        assert(col < getCol());
        const size_t major = MatrixMajor::selectMajor<This>(row, col);
        const size_t minor = MatrixMajor::selectMinor<This>(row, col);
        const size_t from = majorStarts[major];
        const size_t to = majorStarts[major + 1];
        size_t index = 0;
        for (size_t i = from; i < to && index <= minor; ++i) {
            index = minorIndexes[i];
            if (index == minor)
                return elements[i];
        }
        return T(0);
    }

    template<Scalar T, int Major>
    void SparseMatrix<T, Major>::insert(T x, size_t row, size_t col) {
        assert(row < getRow());
        assert(col < getCol());
        assert(!x.isZero() && "[Error]: Zero element is useless");
        const size_t major = MatrixMajor::selectMajor<This>(row, col);
        const size_t minor = MatrixMajor::selectMinor<This>(row, col);
        const size_t from = majorStarts[major];
        const size_t to = majorStarts[major + 1];
        /* Search existing element */ {
            size_t index = 0;
            for (size_t i = from; i < to && index <= minor; ++i) {
                index = minorIndexes[i];
                if (index == minor) {
                    elements[i] = x;
                    return;
                }
            }
        }

        /* If not exist */ {
            size_t insert_to = from;
            while (insert_to < to && minorIndexes[insert_to] < minor)
                ++insert_to;
            elements.insert(insert_to, x);
            minorIndexes.insert(insert_to, minor);
            for (size_t i = major + 1; i < majorStarts.getLength(); ++i)
                ++majorStarts[i];
        }
    }

    template<Scalar T, int Major>
    void SparseMatrix<T, Major>::resize(size_t row, size_t col) {
        zeros();
        majorStarts.resize(MatrixMajor::selectMajor<This>(row, col) + 1, 0);
        maxMinor = MatrixMajor::selectMinor<This>(row, col);
    }

    template<Scalar T, int Major>
    void SparseMatrix<T, Major>::zeros() {
        elements.resize(0);
        minorIndexes.resize(0);
        majorStarts.zeros();
    }

    template<Scalar T, int Major>
    void SparseMatrix<T, Major>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        elements.swap(obj.elements);
        minorIndexes.swap(obj.minorIndexes);
        majorStarts.swap(obj.majorStarts);
        std::swap(maxMinor, obj.maxMinor);
    }

    template<Scalar T, int Major>
    size_t SparseMatrix<T, Major>::getRow() const noexcept {
        return MatrixMajor::isColMatrix<This>() ? getMaxMinor() : getMaxMajor();
    }

    template<Scalar T, int Major>
    size_t SparseMatrix<T, Major>::getCol() const noexcept {
        return MatrixMajor::isColMatrix<This>() ? getMaxMajor() : getMaxMinor();
    }

    template<Scalar T, int Major>
    size_t SparseMatrix<T, Major>::getMaxMajor() const noexcept {
        return majorStarts.getLength() - 1;
    }

    template<Scalar T, int Major>
    size_t SparseMatrix<T, Major>::getMaxMinor() const noexcept {
        return maxMinor;
    }
}

namespace Physica {
    template<Scalar T, int Op>
    class Traits<SparseMatrix<T, Op>> {
    public:
        using ScalarType = T;
    };
}

#include "SparseMatrixImpl/SpMV.h"
