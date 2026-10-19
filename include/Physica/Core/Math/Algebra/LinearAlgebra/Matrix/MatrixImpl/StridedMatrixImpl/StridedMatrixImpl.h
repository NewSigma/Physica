/*
 * Copyright 2026 Weibo He.
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

#include <cstring>
#include "../StridedMatrix.h"

namespace Physica {
    template<class Derived>
    constexpr StridedMatrix<Derived>::StridedMatrix() {
        static_assert(Derived::isRowMatrix() || Derived::isColMatrix(), "[Error]: StridedMatrix requires a known major");
    }

    template<class Derived>
    bool StridedMatrix<Derived>::operator==(this const auto& self, const Matrix auto& other) noexcept {
        if constexpr (isTrivialMinor<decltype(self), decltype(other)>()) {
            if (self.getRow() != other.getRow())
                return false;
            if (self.getCol() != other.getCol())
                return false;
            const size_t maxMinor = self.getMaxMinor();
            for (size_t major = 0; major < self.getMaxMajor(); ++major) {
                const size_t r = self.rowFromMajorMinor(major, 0);
                const size_t c = self.colFromMajorMinor(major, 0);
                if (std::memcmp(self.data_ptr(r, c), other.data_ptr(r, c), maxMinor * sizeof(T)) != 0)
                    return false;
            }
            return true;
        }
        return self.LValueMatrix<Derived>::operator==(other);
    }

    template<class Derived>
    template<ExecutePolicy P>
    void StridedMatrix<Derived>::assign(Matrix auto&& target) const noexcept {
        if constexpr (isTrivialMinor<Derived, std::remove_cvref_t<decltype(target)>>()) {
            const auto& source = Base::getDerived();
            const size_t maxMajor = target.getMaxMajor();
            const size_t maxMinor = target.getMaxMinor();
            for (size_t major = 0; major < maxMajor; ++major) {
                const size_t r = target.rowFromMajorMinor(major, 0);
                const size_t c = target.colFromMajorMinor(major, 0);
                std::memcpy(target.data_ptr(r, c), source.data_ptr(r, c), maxMinor * sizeof(typename StridedMatrix<Derived>::ScalarType));
            }
            return;
        }
        Base::template assign<P>(target);
    }

    template<class Derived>
    constexpr auto StridedMatrix<Derived>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        assert(row < self.getRow());
        assert(col < self.getCol());
        return self.data_handle() + row * self.getRowStride() + col * self.getColStride();
    }

    template<class Derived>
    constexpr auto StridedMatrix<Derived>::row(this auto&& self, size_t r) noexcept {
        using Self = decltype(self);
        return StridedMatrixBlock<Self, 1, Derived::getColAtCompile()>(std::forward<Self>(self), r, 1, 0, self.getCol());
    }

    template<class Derived>
    constexpr auto StridedMatrix<Derived>::col(this auto&& self, size_t c) noexcept {
        using Self = decltype(self);
        return StridedMatrixBlock<Self, Derived::getRowAtCompile(), 1>(std::forward<Self>(self), 0, self.getRow(), c, 1);
    }

    template<class Derived>
    template<size_t Row>
    auto StridedMatrix<Derived>::rows(this auto&& self, size_t fromRow, size_t rowCount) noexcept {
        using Self = decltype(self);
        return StridedMatrixBlock<Self, Row, Derived::getColAtCompile()>(std::forward<Self>(self), fromRow, rowCount, 0, self.getCol());
    }

    template<class Derived>
    template<size_t Row>
    auto StridedMatrix<Derived>::topRows(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return StridedMatrixBlock<Self, Row, Derived::getColAtCompile()>(std::forward<Self>(self), 0, to, 0, self.getCol());
    }

    template<class Derived>
    template<size_t Row>
    auto StridedMatrix<Derived>::bottomRows(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return StridedMatrixBlock<Self, Row, Derived::getColAtCompile()>(std::forward<Self>(self), from, self.getRow() - from, 0, self.getCol());
    }

    template<class Derived>
    template<size_t Col>
    auto StridedMatrix<Derived>::cols(this auto&& self, size_t fromCol, size_t colCount) noexcept {
        using Self = decltype(self);
        return StridedMatrixBlock<Self, Derived::getRowAtCompile(), Col>(std::forward<Self>(self), 0, self.getRow(), fromCol, colCount);
    }

    template<class Derived>
    template<size_t Col>
    auto StridedMatrix<Derived>::leftCols(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return StridedMatrixBlock<Self, Derived::getRowAtCompile(), Col>(std::forward<Self>(self), 0, self.getRow(), 0, to);
    }

    template<class Derived>
    template<size_t Col>
    auto StridedMatrix<Derived>::rightCols(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return StridedMatrixBlock<Self, Derived::getRowAtCompile(), Col>(std::forward<Self>(self), 0, self.getRow(), from, self.getCol() - from);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto StridedMatrix<Derived>::topLeftCorner(this auto&& self, size_t toRow, size_t toCol) noexcept {
        using Self = decltype(self);
        return StridedMatrixBlock<Self, Row, Col>(std::forward<Self>(self), 0, toRow, 0, toCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto StridedMatrix<Derived>::topLeftCorner(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return StridedMatrixBlock<Self, Row, Col>(std::forward<Self>(self), 0, to, 0, to);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto StridedMatrix<Derived>::topRightCorner(this auto&& self, size_t toRow, size_t fromCol) noexcept {
        using Self = decltype(self);
        return StridedMatrixBlock<Self, Row, Col>(std::forward<Self>(self), 0, toRow, fromCol, self.getRow() - fromCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto StridedMatrix<Derived>::bottomLeftCorner(this auto&& self, size_t fromRow, size_t toCol) noexcept {
        using Self = decltype(self);
        return StridedMatrixBlock<Self, Row, Col>(std::forward<Self>(self), fromRow, self.getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto StridedMatrix<Derived>::bottomRightCorner(this auto&& self, size_t fromRow, size_t fromCol) noexcept {
        using Self = decltype(self);
        return StridedMatrixBlock<Self, Row, Col>(std::forward<Self>(self), fromRow, self.getRow() - fromRow, fromCol, self.getCol() - fromCol);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto StridedMatrix<Derived>::bottomRightCorner(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return StridedMatrixBlock<Self, Row, Col>(std::forward<Self>(self), from, self.getRow() - from, from, self.getCol() - from);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto StridedMatrix<Derived>::block(this auto&& self, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        using Self = decltype(self);
        return StridedMatrixBlock<Self, Row, Col>(std::forward<Self>(self), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    void StridedMatrix<Derived>::zeros() noexcept {
        if constexpr (isTrivialMinor<Derived, Derived>()) {
            auto& self = Base::getDerived();
            const size_t maxMajor = self.getMaxMajor();
            const size_t maxMinor = self.getMaxMinor();
            for (size_t major = 0; major < maxMajor; ++major) {
                const size_t r = self.rowFromMajorMinor(major, 0);
                const size_t c = self.colFromMajorMinor(major, 0);
                std::memset(self.data_ptr(r, c), 0, maxMinor * sizeof(typename Base::ScalarType));
            }
            return;
        }
        Base::zeros();
    }

    template<class Derived>
    constexpr size_t StridedMatrix<Derived>::getRowStride() const noexcept {
        if constexpr (Derived::getRowStrideAtCompile() != Dynamic)
            return Derived::getRowStrideAtCompile();
        else
            return Base::getDerived().getRowStride();
    }

    template<class Derived>
    constexpr size_t StridedMatrix<Derived>::getColStride() const noexcept {
        if constexpr (Derived::getColStrideAtCompile() != Dynamic)
            return Derived::getColStrideAtCompile();
        else
            return Base::getDerived().getColStride();
    }

    template<class Derived>
    constexpr size_t StridedMatrix<Derived>::getMajorStride() const noexcept {
        if constexpr (Derived::isRowMatrix())
            return getRowStride();
        else
            return getColStride();
    }

    template<class Derived>
    constexpr size_t StridedMatrix<Derived>::getMinorStride() const noexcept {
        if constexpr (Derived::isRowMatrix())
            return getColStride();
        else
            return getRowStride();
    }

    template<class Derived>
    auto StridedMatrix<Derived>::data_handle() noexcept {
        return Base::getDerived().data_handle();
    }

    template<class Derived>
    auto StridedMatrix<Derived>::data_handle() const noexcept {
        return Base::getDerived().data_handle();
    }

    template<class Derived>
    __host__ __device__ consteval size_t StridedMatrix<Derived>::getMajorStrideAtCompile() noexcept {
        return Derived::isRowMatrix() ? Derived::getRowStrideAtCompile() : Derived::getColStrideAtCompile();
    }

    template<class Derived>
    __host__ __device__ consteval size_t StridedMatrix<Derived>::getMinorStrideAtCompile() noexcept {
        return Derived::isRowMatrix() ? Derived::getColStrideAtCompile() : Derived::getRowStrideAtCompile();
    }

    template<class Derived>
    template<class M1, class M2>
    consteval bool StridedMatrix<Derived>::isTrivialMinor() noexcept {
        using U1 = std::remove_cvref_t<M1>;
        using U2 = std::remove_cvref_t<M2>;
        using T1 = U1::ScalarType;
        using T2 = U2::ScalarType;
        if constexpr (U1::isStrided() && U2::isStrided()) {
            return U1::getMinorStrideAtCompile() == 1 && U2::getMinorStrideAtCompile() == 1
                && MatrixMajor::isSameMajor<M1, M2>()
                && !U1::isBothMajor() && !U2::isBothMajor()
                && !U1::isDiffable() && !U2::isDiffable()
                && std::same_as<T1, T2>
                && std::is_trivially_copyable_v<T1>;
        }
        else
            return false;
    }
}
