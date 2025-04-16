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

#include <cassert>
#include <forward_list>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RValueMatrix.h"
#include "FormatedMatrix.h"

namespace Physica {
    template<class Derived>
    template<Matrix M>
    void RValueMatrix<Derived>::assign(M& target) const {
        assert(getRow() == target.getRow() && "[Error]: Dimensions do not match");
        assert(getCol() == target.getCol() && "[Error]: Dimensions do not match");
        assign_check(target);

        const size_t maxMajor = target.getMaxMajor();
        const size_t maxMinor = target.getMaxMinor();
        for (size_t i = 0; i < maxMajor; ++i)
            for (size_t j = 0; j < maxMinor; ++j)
                target.refFromMajorMinor(i, j) = calc(target.rowFromMajorMinor(i, j), target.colFromMajorMinor(i, j));
    }

    template<class Derived>
    template<Matrix M>
    void RValueMatrix<Derived>::assign_add(M& target) const {
        assert(getRow() == target.getRow() && "[Error]: Dimensions do not match");
        assert(getCol() == target.getCol() && "[Error]: Dimensions do not match");
        assign_check(target);

        const size_t maxMajor = target.getMaxMajor();
        const size_t maxMinor = target.getMaxMinor();
        for (size_t i = 0; i < maxMajor; ++i)
            for (size_t j = 0; j < maxMinor; ++j)
                target.refFromMajorMinor(i, j) += calc(target.rowFromMajorMinor(i, j), target.colFromMajorMinor(i, j));
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::row(size_t r) noexcept {
        return RowVector(Base::getDerived(), r, 0, getCol());
    }

    template<class Derived>
    inline const auto RValueMatrix<Derived>::row(size_t r) const noexcept {
        return RowVector(Base::getConstCastDerived(), r, 0, getCol());
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::col(size_t c) noexcept {
        return ColVector(Base::getDerived(), 0, getRow(), c);
    }

    template<class Derived>
    inline const auto RValueMatrix<Derived>::col(size_t c) const noexcept {
        return ColVector(Base::getConstCastDerived(), 0, getRow(), c);
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::rows(size_t fromRow, size_t rowCount) noexcept {
        return BlockType(Base::getDerived(), fromRow, rowCount, 0, getCol());
    }

    template<class Derived>
    inline const auto RValueMatrix<Derived>::rows(size_t fromRow, size_t rowCount) const noexcept {
        return BlockType(Base::getConstCastDerived(), fromRow, rowCount, 0, getCol());
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::topRows(size_t to) noexcept {
        return BlockType(Base::getDerived(), 0, to, 0, getCol());
    }

    template<class Derived>
    inline const auto RValueMatrix<Derived>::topRows(size_t to) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, to, 0, getCol());
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::bottomRows(size_t from) noexcept {
        return BlockType(Base::getDerived(), from, getRow() - from, 0, getCol());
    }

    template<class Derived>
    inline const auto RValueMatrix<Derived>::bottomRows(size_t from) const noexcept {
        return BlockType(Base::getConstCastDerived(), from, getRow() - from, 0, getCol());
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::cols(size_t fromCol, size_t colCount) noexcept {
        return BlockType(Base::getDerived(), 0, getRow(), fromCol, colCount);
    }

    template<class Derived>
    inline const auto RValueMatrix<Derived>::cols(size_t fromCol, size_t colCount) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, getRow(), fromCol, colCount);
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::leftCols(size_t to) noexcept {
        return BlockType(Base::getDerived(), 0, getRow(), 0, to);
    }

    template<class Derived>
    inline const auto RValueMatrix<Derived>::leftCols(size_t to) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, getRow(), 0, to);
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::rightCols(size_t from) noexcept {
        return BlockType(Base::getDerived(), 0, getRow(), from, getCol() - from);
    }

    template<class Derived>
    inline const auto RValueMatrix<Derived>::rightCols(size_t from) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, getRow(), from, getCol() - from);
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::topLeftCorner(size_t toRow, size_t toCol) noexcept {
        return BlockType(Base::getDerived(), 0, toRow, 0, toCol);
    }

    template<class Derived>
    inline const auto RValueMatrix<Derived>::topLeftCorner(size_t toRow, size_t toCol) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, toRow, 0, toCol);
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::topLeftCorner(size_t to) noexcept {
        return BlockType(Base::getDerived(), 0, to, 0, to);
    }

    template<class Derived>
    inline const auto RValueMatrix<Derived>::topLeftCorner(size_t to) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, to, 0, to);
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::topRightCorner(size_t toRow, size_t fromCol) noexcept {
        return BlockType(Base::getDerived(), 0, toRow, fromCol, getRow() - fromCol);
    }

    template<class Derived>
    inline const auto RValueMatrix<Derived>::topRightCorner(size_t toRow, size_t fromCol) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, toRow, fromCol, getRow() - fromCol);
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::bottomLeftCorner(size_t fromRow, size_t toCol) noexcept {
        return BlockType(Base::getDerived(), fromRow, getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    inline const auto RValueMatrix<Derived>::bottomLeftCorner(size_t fromRow, size_t toCol) const noexcept {
        return BlockType(Base::getConstCastDerived(), fromRow, getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::bottomRightCorner(size_t fromRow, size_t fromCol) noexcept {
        return BlockType(Base::getDerived(), fromRow, getRow() - fromRow, fromCol, getCol() - fromCol);
    }

    template<class Derived>
    inline const auto RValueMatrix<Derived>::bottomRightCorner(size_t fromRow, size_t fromCol) const noexcept {
        return BlockType(Base::getConstCastDerived(), fromRow, getRow() - fromRow, fromCol, getCol() - fromCol);
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::bottomRightCorner(size_t from) noexcept {
        return BlockType(Base::getDerived(), from, getRow() - from, from, getCol() - from);
    }

    template<class Derived>
    inline const auto RValueMatrix<Derived>::bottomRightCorner(size_t from) const noexcept {
        return BlockType(Base::getConstCastDerived(), from, getRow() - from, from, getCol() - from);
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        return BlockType(Base::getDerived(), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    inline const auto RValueMatrix<Derived>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const noexcept {
        return BlockType(Base::getConstCastDerived(), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::diag() noexcept {
        return DiagVector<Derived, false>(Base::getDerived());
    }

    template<class Derived>
    inline const auto RValueMatrix<Derived>::diag() const noexcept {
        return DiagVector<Derived, false>(Base::getConstCastDerived());
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::calcFromMajorMinor(size_t major, size_t minor) const {
        return calc(rowFromMajorMinor(major, minor), colFromMajorMinor(major, minor));
    }

    template<class Derived>
    template<Matrix M1, Matrix M2>
    void RValueMatrix<Derived>::reverse(const M1&, const M2& grad) const noexcept requires(isReverseDiff) {
        Base::getDerived().reverse(grad);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::max() const -> T {
        T result;
        if constexpr (MatrixOption::isColMatrix<This>()) {
            result = col(0).max();
            for (size_t i = 1; i < getCol(); ++i) {
                T temp = col(i).max();
                if (temp > result)
                    result = temp;
            }
        }
        else {
            result = row(0).max();
            for (size_t i = 1; i < getRow(); ++i) {
                T temp = row(i).max();
                if (temp > result)
                    result = temp;
            }
        }
        return result;
    }

    template<class Derived>
    auto RValueMatrix<Derived>::min() const -> T {
        T result;
        if constexpr (MatrixOption::isColMatrix<This>()) {
            result = col(0).min();
            for (size_t i = 1; i < getCol(); ++i) {
                T temp = col(i).min();
                if (temp < result)
                    result = temp;
            }
        }
        else {
            result = row(0).min();
            for (size_t i = 1; i < getRow(); ++i) {
                T temp = row(i).min();
                if (temp < result)
                    result = temp;
            }
        }
        return result;
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::sum() const -> CoDiff<T> {
        if constexpr (isReverseDiff) {
            Tv v = 0;
            std::forward_list<CoDiff<T>> elems{};
            for (size_t major = 0; major < getMaxMajor(); ++major) {
                size_t minor = 0;
                auto list = co_for([&]{ return minor < getMaxMinor(); }, [&]{ ++minor; }, [&]{
                    auto elem = calcFromMajorMinor(major, minor);
                    v += elem.value();
                    return elem;
                });
                elems.merge(std::move(list));
            }
            auto result = co_yield std::move(v);
            for (auto& elem : elems)
                elem.reverse(result.grad());
        }
        else {
            T result = 0;
            for (size_t major = 0; major < getMaxMajor(); ++major)
                for (size_t minor = 0; minor < getMaxMinor(); ++minor)
                    result += calcFromMajorMinor(major, minor);
            co_return std::move(result);
        }
    }

    template<class Derived>
    auto RValueMatrix<Derived>::sum_cols() const {
        return MatrixSum<Derived, true>(Base::getDerived());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::mean() const -> T {
        return sum() / Trv(getRow() * getCol());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::trace() const -> T {
        assert(getRow() == getCol());
        T result = T(0);
        for (size_t i = 0; i < getRow(); ++i)
            result += calc(i, i);
        return result;
    }

    template<class Derived>
    auto RValueMatrix<Derived>::lnSumExp() const -> CoDiff<T> {
        return flatten().lnSumExp();
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::format() const noexcept {
        return FormatedMatrix<Derived>(Base::getDerived());
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::transpose() const noexcept {
        return Transpose<Derived>(Base::getDerived());
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::conjugate() const noexcept {
        return Conjugate<Derived>(Base::getDerived());
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::hermite() const noexcept {
        return Hermite<Derived>(Base::getDerived());
    }

    template<class Derived>
    inline auto RValueMatrix<Derived>::flatten() const noexcept {
        return FlattenR<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::reals() const noexcept -> RealsRtnTy {
        return RealsRtnTy(Base::getDerived());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::imags() const noexcept {
        return ImagMatrix<Derived>(Base::getDerived());
    }
    
    template<class Derived>
    auto RValueMatrix<Derived>::squaredNorms() const noexcept {
        return SquaredNormMatrix<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::norms() const noexcept {
        return NormMatrix<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::values() const noexcept -> ValuesRtnTy {
        return ValuesRtnTy(Base::getDerived());
    }

    template<class Derived>
    template<int GradOrder>
    auto RValueMatrix<Derived>::grads() const noexcept {
        if constexpr (isReverseDiff)
            return Base::getDerived().template grads<GradOrder>();
        else
            return grads_impl<GradOrder>();
    }

    template<class Derived>
    template<int GradOrder>
    auto RValueMatrix<Derived>::grads_impl() const noexcept {
        return GradMatrix<Derived, GradOrder>(Base::getDerived());
    }

    template<class Derived>
    bool RValueMatrix<Derived>::isSymm() const noexcept {
        if (getRow() != getCol())
            return false;
        const size_t order = getRow();
        for (size_t r = 0; r < order; ++r)
            for (size_t c = r + 1; c < order; ++c)
                if (calc(r, c) != calc(c, r))
                    return false;
        return true;
    }

    template<class Derived>
    bool RValueMatrix<Derived>::isFinite() const noexcept {
        for (size_t r = 0; r < getRow(); ++r)
            for (size_t c = 0; c < getCol(); ++c)
                if (!calc(r, c).isFinite())
                    return false;
        return true;
    }

    template<class Derived>
    template<Matrix M>
    __host__ __device__ void RValueMatrix<Derived>::assign_check(const M&) noexcept {
        static_assert(RowAtCompile == M::RowAtCompile || RowAtCompile == Dynamic || M::RowAtCompile == Dynamic, "[Error]: Row mismatch between two matrix");
        static_assert(ColAtCompile == M::ColAtCompile || ColAtCompile == Dynamic || M::ColAtCompile == Dynamic, "[Error]: Col mismatch between two matrix");
        static_assert(!isComplex || M::isComplex, "[Error]: Assign a complex matrix to real matrix discards imag part");
    }

    template<Matrix T1, Matrix T2>
    bool matrixNear(const T1& m1, const T2& m2, double precision) {
        assert(m1.getRow() == m2.getRow());
        assert(m1.getCol() == m2.getCol());
        for (size_t i = 0; i < m1.getCol(); ++i)
            for (size_t j = 0; j < m1.getRow(); ++j)
                if (!scalarNear(m1.calc(j, i), m2.calc(j, i), precision))
                    return false;
        return true;
    }

    template<Matrix T>
    bool operator==(const T& m1, const T& m2) {
        if (m1.getRow() != m2.getRow())
            return false;
        if (m1.getCol() != m2.getCol())
            return false;
        for (size_t major = 0; major < m1.getMaxMajor(); ++major)
            for (size_t minor = 0; minor < m1.getMaxMinor(); ++minor)
                if (m1.calcFromMajorMinor(major, minor) != m2.calcFromMajorMinor(major, minor))
                    return false;
        return true;
    }
}
