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
    template<ExecutePolicy P>
    void RValueMatrix<Derived>::assign(Matrix auto& target) const {
        assert(getRow() == target.getRow() && "[Error]: Dimensions do not match");
        assert(getCol() == target.getCol() && "[Error]: Dimensions do not match");
        assert_assign(target);

        const size_t maxMajor = target.getMaxMajor();
        const size_t maxMinor = target.getMaxMinor();
        for (size_t i = 0; i < maxMajor; ++i)
            for (size_t j = 0; j < maxMinor; ++j)
                target.refFromMajorMinor(i, j) = calc(target.rowFromMajorMinor(i, j), target.colFromMajorMinor(i, j));
    }

    template<class Derived>
    void RValueMatrix<Derived>::assign_add(Matrix auto& target) const {
        assert(getRow() == target.getRow() && "[Error]: Dimensions do not match");
        assert(getCol() == target.getCol() && "[Error]: Dimensions do not match");
        assert_assign(target);

        const size_t maxMajor = target.getMaxMajor();
        const size_t maxMinor = target.getMaxMinor();
        for (size_t i = 0; i < maxMajor; ++i)
            for (size_t j = 0; j < maxMinor; ++j)
                target.refFromMajorMinor(i, j) += calc(target.rowFromMajorMinor(i, j), target.colFromMajorMinor(i, j));
    }

    template<class Derived>
    void RValueMatrix<Derived>::assert_assign(const Matrix auto& target) const noexcept {
        static_assert_assign(target);
        assert(getRow() == target.getRow() && "[Error]: Dimensions do not match");
        assert(getCol() == target.getCol() && "[Error]: Dimensions do not match");
    }

    template<class Derived>
    decltype(auto) RValueMatrix<Derived>::calcFromMajorMinor(size_t major, size_t minor) const {
        return calc(rowFromMajorMinor(major, minor), colFromMajorMinor(major, minor));
    }

    template<class Derived>
    void RValueMatrix<Derived>::reverse(const Matrix auto&, const Matrix auto& grad) const noexcept requires(isReverseDiff) {
        Base::getDerived().reverse(grad);
    }

    template<class Derived>
    void RValueMatrix<Derived>::resize(const Matrix auto& m, auto&&... args) {
        resize(m.getRow(), m.getCol(), std::forward<decltype(args)>(args)...);
    }

    template<class Derived>
    decltype(auto) RValueMatrix<Derived>::resize(size_t r, size_t c, auto&&... args) {
        return Base::getDerived().resize(r, c, std::forward<decltype(args)>(args)...);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::row(size_t r) noexcept {
        return RowVector(Base::getDerived(), r, 0, getCol());
    }

    template<class Derived>
    const auto RValueMatrix<Derived>::row(size_t r) const noexcept {
        return RowVector(Base::getConstCastDerived(), r, 0, getCol());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::col(size_t c) noexcept {
        return ColVector(Base::getDerived(), 0, getRow(), c);
    }

    template<class Derived>
    const auto RValueMatrix<Derived>::col(size_t c) const noexcept {
        return ColVector(Base::getConstCastDerived(), 0, getRow(), c);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::rows(size_t fromRow, size_t rowCount) noexcept {
        return BlockType(Base::getDerived(), fromRow, rowCount, 0, getCol());
    }

    template<class Derived>
    const auto RValueMatrix<Derived>::rows(size_t fromRow, size_t rowCount) const noexcept {
        return BlockType(Base::getConstCastDerived(), fromRow, rowCount, 0, getCol());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::topRows(size_t to) noexcept {
        return BlockType(Base::getDerived(), 0, to, 0, getCol());
    }

    template<class Derived>
    const auto RValueMatrix<Derived>::topRows(size_t to) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, to, 0, getCol());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::bottomRows(size_t from) noexcept {
        return BlockType(Base::getDerived(), from, getRow() - from, 0, getCol());
    }

    template<class Derived>
    const auto RValueMatrix<Derived>::bottomRows(size_t from) const noexcept {
        return BlockType(Base::getConstCastDerived(), from, getRow() - from, 0, getCol());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::cols(size_t fromCol, size_t colCount) noexcept {
        return BlockType(Base::getDerived(), 0, getRow(), fromCol, colCount);
    }

    template<class Derived>
    const auto RValueMatrix<Derived>::cols(size_t fromCol, size_t colCount) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, getRow(), fromCol, colCount);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::leftCols(size_t to) noexcept {
        return BlockType(Base::getDerived(), 0, getRow(), 0, to);
    }

    template<class Derived>
    const auto RValueMatrix<Derived>::leftCols(size_t to) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, getRow(), 0, to);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::rightCols(size_t from) noexcept {
        return BlockType(Base::getDerived(), 0, getRow(), from, getCol() - from);
    }

    template<class Derived>
    const auto RValueMatrix<Derived>::rightCols(size_t from) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, getRow(), from, getCol() - from);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::topLeftCorner(size_t toRow, size_t toCol) noexcept {
        return BlockType(Base::getDerived(), 0, toRow, 0, toCol);
    }

    template<class Derived>
    const auto RValueMatrix<Derived>::topLeftCorner(size_t toRow, size_t toCol) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, toRow, 0, toCol);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::topLeftCorner(size_t to) noexcept {
        return BlockType(Base::getDerived(), 0, to, 0, to);
    }

    template<class Derived>
    const auto RValueMatrix<Derived>::topLeftCorner(size_t to) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, to, 0, to);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::topRightCorner(size_t toRow, size_t fromCol) noexcept {
        return BlockType(Base::getDerived(), 0, toRow, fromCol, getRow() - fromCol);
    }

    template<class Derived>
    const auto RValueMatrix<Derived>::topRightCorner(size_t toRow, size_t fromCol) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, toRow, fromCol, getRow() - fromCol);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::bottomLeftCorner(size_t fromRow, size_t toCol) noexcept {
        return BlockType(Base::getDerived(), fromRow, getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    const auto RValueMatrix<Derived>::bottomLeftCorner(size_t fromRow, size_t toCol) const noexcept {
        return BlockType(Base::getConstCastDerived(), fromRow, getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::bottomRightCorner(size_t fromRow, size_t fromCol) noexcept {
        return BlockType(Base::getDerived(), fromRow, getRow() - fromRow, fromCol, getCol() - fromCol);
    }

    template<class Derived>
    const auto RValueMatrix<Derived>::bottomRightCorner(size_t fromRow, size_t fromCol) const noexcept {
        return BlockType(Base::getConstCastDerived(), fromRow, getRow() - fromRow, fromCol, getCol() - fromCol);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::bottomRightCorner(size_t from) noexcept {
        return BlockType(Base::getDerived(), from, getRow() - from, from, getCol() - from);
    }

    template<class Derived>
    const auto RValueMatrix<Derived>::bottomRightCorner(size_t from) const noexcept {
        return BlockType(Base::getConstCastDerived(), from, getRow() - from, from, getCol() - from);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        return BlockType(Base::getDerived(), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    const auto RValueMatrix<Derived>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const noexcept {
        return BlockType(Base::getConstCastDerived(), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::diag() noexcept {
        return DiagVector<Derived, false>(Base::getDerived());
    }

    template<class Derived>
    const auto RValueMatrix<Derived>::diag() const noexcept {
        return DiagVector<Derived, false>(Base::getConstCastDerived());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::triu() noexcept {
        return TrigUpper<Derived>(Base::getDerived());
    }

    template<class Derived>
    const auto RValueMatrix<Derived>::triu() const noexcept {
        return TrigUpper<Derived>(Base::getConstCastDerived());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::tril() noexcept {
        return TrigLower<Derived>(Base::getDerived());
    }

    template<class Derived>
    const auto RValueMatrix<Derived>::tril() const noexcept {
        return TrigLower<Derived>(Base::getConstCastDerived());
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
    auto RValueMatrix<Derived>::sum() const -> CoDiff<T> {
        if constexpr (isReverseDiff) {
            auto result = co_yield values().sum();
            Base::getDerived().reverse(result.grad());
        }
        else {
            T result = 0;
            for (size_t major = 0; major < getMaxMajor(); ++major)
                result += MatrixOption::isColMatrix<Derived>() ? col(major).sum() : row(major).sum();
            co_return std::move(result);
        }
    }

    template<class Derived>
    auto RValueMatrix<Derived>::sum_rows() const {
        return MatrixSum<Derived, false>(Base::getDerived());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::sum_cols() const {
        return MatrixSum<Derived, true>(Base::getDerived());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::mean() const -> CoDiff<T> {
        return sum() / Trv(getSize());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::trace() const -> T {
        assert(isSquare());
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
    auto RValueMatrix<Derived>::det() const -> CoDiff<T> {
        assert(isSquare() && "[Error]: Determinate requires square matrix");
        constexpr size_t Order = RowAtCompile > ColAtCompile ? RowAtCompile : ColAtCompile;
        if constexpr (Order == 1)
            return calc(0, 0);
        else if constexpr (Order == 2)
            return calc(0, 0) * calc(1, 1) - calc(0, 1) * calc(1, 0);
        else if constexpr (Order == 3)
            return calc(0, 0) * (calc(1, 1) * calc(2, 2) - calc(1, 2) * calc(2, 1))
                 + calc(0, 1) * (calc(1, 2) * calc(2, 0) - calc(1, 0) * calc(2, 2))
                 + calc(0, 2) * (calc(1, 0) * calc(2, 1) - calc(1, 1) * calc(2, 0));
        else
            return QRDecomp<T, false>(Base::getDerived()).det();
    }

    template<class Derived>
    auto RValueMatrix<Derived>::lnAbsDet() const -> Tr {
        QRDecomp<T, false> qr(Base::getDerived());
        return ln(abs(qr.getMatrixR().diag())).sum();
    }

    template<class Derived>
    auto RValueMatrix<Derived>::sgndet() const {
        static_assert(!T::isComplex, "[Error]: sgndet() is not well defined for complex matrix");
        return Base::getDerived().sgndet();
    }

    template<class Derived>
    auto RValueMatrix<Derived>::format() const noexcept {
        return FormatedMatrix<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::inv() const noexcept {
        return Inverse<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::pinv() const noexcept {
        return PseudoInverse<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::transpose() const noexcept {
        return Transpose<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::conjugate() const noexcept -> ConjugateRtnTy {
        if constexpr (isComplex)
            return Conjugate<Derived>(Base::getDerived());
        else
            return Base::getDerived();
    }

    template<class Derived>
    auto RValueMatrix<Derived>::hermite() const noexcept {
        return Hermite<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::flatten() const noexcept {
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
    bool RValueMatrix<Derived>::isOverdetermined() const noexcept {
        return getRow() > getCol();
    }

    template<class Derived>
    bool RValueMatrix<Derived>::isUnderdetermined() const noexcept {
        return getRow() < getCol();
    }

    template<class Derived>
    bool RValueMatrix<Derived>::isSquare() const noexcept {
        return getRow() == getCol();
    }

    template<class Derived>
    bool RValueMatrix<Derived>::isSymm() const noexcept {
        if constexpr (MatrixOption::isSymmMatrix<Derived>())
            return true;
        if (!isSquare())
            return false;

        const size_t order = getRow();
        for (size_t r = 0; r < order; ++r)
            for (size_t c = r + 1; c < order; ++c)
                if (calc(r, c) != calc(c, r))
                    return false;
        return true;
    }

    template<class Derived>
    bool RValueMatrix<Derived>::isHermite() const noexcept {
        if constexpr (MatrixOption::isHermiteMatrix<Derived>())
            return true;
        if (!isSquare())
            return false;

        const size_t order = getRow();
        for (size_t r = 0; r < order; ++r)
            for (size_t c = r + 1; c < order; ++c)
                if (calc(r, c) != calc(c, r).conjugate())
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
    __host__ __device__ void RValueMatrix<Derived>::static_assert_assign(const Matrix auto& target) noexcept {
        using M = std::remove_cvref_t<decltype(target)>;
        static_assert(RowAtCompile == M::RowAtCompile || RowAtCompile == Dynamic || M::RowAtCompile == Dynamic, "[Error]: Row mismatch between two matrix");
        static_assert(ColAtCompile == M::ColAtCompile || ColAtCompile == Dynamic || M::ColAtCompile == Dynamic, "[Error]: Col mismatch between two matrix");
        static_assert(!isComplex || M::isComplex, "[Error]: Assign a complex matrix to real matrix discards imag part");
        static_assert(!isDiffable || M::isDiffable, "[Error]: Assign a diffable matrix to normal matrix discards grads");
    }
}

#include "Sum.h"
#include "Inverse.h"
#include "PseudoInverse.h"
#include "Transpose.h"
#include "Conjugate.h"
#include "Hermite.h"
#include "Flatten.h"
#include "Trig/Trig.h"
#include "MatrixConvert.h"
#include "ReshapedVector.h"
