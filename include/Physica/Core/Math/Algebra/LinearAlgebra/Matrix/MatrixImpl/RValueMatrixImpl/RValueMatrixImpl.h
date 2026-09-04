/*
 * Copyright 2021-2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RValueMatrix.h"
#include "FormatedMatrix.h"

namespace Physica {
    template<class Derived, Scalar ScalarT>
    bool RValueMatrix<Derived, ScalarT>::operator!=(this const auto& self, const Matrix auto& other) noexcept {
        return !(self == other);
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::operator*(this auto&& self, Vector auto&& v) noexcept {
        assert(self.getCol() == v.getLength());
        using Self = decltype(self);
        using V = decltype(v);
        static_assert(!is_device_obj<V>::value, "[Error]: host-device mismatch");
        if constexpr (self.getRowAtCompile() == 1)
            return std::forward<Self>(self).row(0) * v;
        else
            return GEMV<Self, V&&>(std::forward<Self>(self), std::forward<V>(v));
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::operator*(this auto&& self, Matrix auto&& m) noexcept {
        using Self = decltype(self);
        using M = decltype(m);
        static_assert(!is_device_obj<M>::value, "[Error]: host-device mismatch");
        constexpr bool ColVectorLHS = self.getColAtCompile() == 1;
        constexpr bool RowVectorLHS = self.getRowAtCompile() == 1;
        constexpr bool ColVectorRHS = m.getColAtCompile() == 1;
        constexpr bool RowVectorRHS = m.getRowAtCompile() == 1;
        if constexpr (RowVectorLHS)
            return (std::forward<M>(m).transpose() * std::forward<Self>(self).row(0)).transpose();
        else if constexpr (ColVectorRHS)
            return std::forward<Self>(self) * std::forward<M>(m).col(0);
        else if constexpr (ColVectorLHS || RowVectorRHS) {
            assert(self.getCol() == m.getRow());
            return std::forward<Self>(self).col(0) * std::forward<M>(m);
        }
        else
            return GEMM<Self, M>(std::forward<Self>(self), std::forward<M>(m));
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::operator-(this auto&& self) noexcept {
        using Self = decltype(self);
        return MatrixExpr<ExprID::Minus, Self>(std::forward<Self>(self));
    }

    template<class Derived, Scalar ScalarT>
    template<ExecutePolicy P>
    void RValueMatrix<Derived, ScalarT>::assign(Matrix auto&& target) const noexcept {
        if constexpr (!isDiffable() && target.isDiffable()) {
            Base::getDerived().assign(target.values());
            target.zero_grad();
        }
        else {
            target.assert_assign(Base::getDerived());

            const size_t maxMajor = target.getMaxMajor();
            const size_t maxMinor = target.getMaxMinor();
            for (size_t i = 0; i < maxMajor; ++i)
                for (size_t j = 0; j < maxMinor; ++j)
                    target.refFromMajorMinor(i, j) = calc(target.rowFromMajorMinor(i, j), target.colFromMajorMinor(i, j));
        }
    }

    template<class Derived, Scalar ScalarT>
    template<ExecutePolicy P>
    void RValueMatrix<Derived, ScalarT>::assign_base(Matrix auto&& target) const noexcept {
        assign<P>(target);
    }

    template<class Derived, Scalar ScalarT>
    template<ExecutePolicy P>
    void RValueMatrix<Derived, ScalarT>::assign_add(Matrix auto&& target) const noexcept {
        if constexpr (!isDiffable() && target.isDiffable())
            Base::getDerived().template assign_add<P>(target.values());
        else {
            target.assert_assign(Base::getDerived());

            const size_t maxMajor = target.getMaxMajor();
            const size_t maxMinor = target.getMaxMinor();
            for (size_t i = 0; i < maxMajor; ++i)
                for (size_t j = 0; j < maxMinor; ++j)
                    target.refFromMajorMinor(i, j) += calc(target.rowFromMajorMinor(i, j), target.colFromMajorMinor(i, j));
        }
    }

    template<class Derived, Scalar ScalarT>
    void RValueMatrix<Derived, ScalarT>::assert_assign(const Matrix auto& source) const noexcept {
        static_assert_assign(source);
        if constexpr (std::is_same<Derived, std::remove_cvref_t<decltype(source)>>::value)
            assert(this != &source && "[Error]: Self assign is likely a bug");

        constexpr size_t Row1 = Derived::getRowAtCompile();
        constexpr size_t Row2 = source.getRowAtCompile();
        if constexpr (Row1 == Dynamic || Row2 == Dynamic)
            assert(getRow() == source.getRow() && "[Error]: Dimensions do not match");

        constexpr size_t Col1 = Derived::getColAtCompile();
        constexpr size_t Col2 = source.getColAtCompile();
        if constexpr (Col1 == Dynamic || Col2 == Dynamic)
            assert(getCol() == source.getCol() && "[Error]: Dimensions do not match");

        constexpr size_t Size1 = Derived::getSizeAtCompile();
        constexpr size_t Size2 = source.getSizeAtCompile();
        if constexpr (Size1 == Dynamic || Size2 == Dynamic)
            assert(getSize() > 0 && "[Error]: Assign a empty matrix is not allowed");
    }

    template<class Derived, Scalar ScalarT>
    void RValueMatrix<Derived, ScalarT>::assert_assign_lapack(const Matrix auto& source) const noexcept {
        static_assert(Internal::EnableLAPACK<Derived, decltype(source)>::value, "[Error]: Invalid expr for LAPACK");
        assert_assign(source);
    }

    template<class Derived, Scalar ScalarT>
    decltype(auto) RValueMatrix<Derived, ScalarT>::calc(size_t row, size_t col) const {
        return Base::getDerived().calc(row, col);
    }

    template<class Derived, Scalar ScalarT>
    decltype(auto) RValueMatrix<Derived, ScalarT>::calc_value(size_t row, size_t col) const {
        return Base::getDerived().values().calc(row, col);
    }

    template<class Derived, Scalar ScalarT>
    decltype(auto) RValueMatrix<Derived, ScalarT>::calcFromMajorMinor(size_t major, size_t minor) const {
        return calc(rowFromMajorMinor(major, minor), colFromMajorMinor(major, minor));
    }

    template<class Derived, Scalar ScalarT>
    void RValueMatrix<Derived, ScalarT>::reverse(this const auto& self, const Matrix auto& grad) noexcept {
        static_assert(isReverseDiff());
        self.grads().assert_assign(grad);
        for (size_t major = 0; major < self.getMaxMajor(); ++major) {
            for (size_t minor = 0; minor < self.getMaxMinor(); ++minor) {
                size_t r = rowFromMajorMinor(major, minor);
                size_t c = colFromMajorMinor(major, minor);
                self.calc(major, minor).reverse(grad.calc(r, c));
            }
        }
    }

    template<class Derived, Scalar ScalarT>
    void RValueMatrix<Derived, ScalarT>::reverse(this const auto& self, const Matrix auto&, const Matrix auto& grad) noexcept {
        static_assert(isReverseDiff());
        self.reverse(grad);
    }

    template<class Derived, Scalar ScalarT>
    void RValueMatrix<Derived, ScalarT>::resize(this auto& self, size_t order) {
        self.resize(order, order);
    }

    template<class Derived, Scalar ScalarT>
    void RValueMatrix<Derived, ScalarT>::resize(this auto& self, const Matrix auto& m, auto&&... args) {
        self.resize(m.getRow(), m.getCol(), std::forward<decltype(args)>(args)...);
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::resize(this auto& self, size_t r, size_t c, auto&&... args) {
        self.resize(r, c, std::forward<decltype(args)>(args)...);
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::row(this auto&& self, size_t r) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self, 1, Dynamic>(std::forward<Self>(self), r, 1, 0, self.getCol());
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::col(this auto&& self, size_t c) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self, Dynamic, 1>(std::forward<Self>(self), 0, self.getRow(), c, 1);
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::rows(this auto&& self, size_t fromRow, size_t rowCount) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), fromRow, rowCount, 0, self.getCol());
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::topRows(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), 0, to, 0, self.getCol());
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::bottomRows(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), from, self.getRow() - from, 0, self.getCol());
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::cols(this auto&& self, size_t fromCol, size_t colCount) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), 0, self.getRow(), fromCol, colCount);
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::leftCols(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), 0, self.getRow(), 0, to);
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::rightCols(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), 0, self.getRow(), from, self.getCol() - from);
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::topLeftCorner(this auto&& self, size_t toRow, size_t toCol) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), 0, toRow, 0, toCol);
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::topLeftCorner(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), 0, to, 0, to);
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::topRightCorner(this auto&& self, size_t toRow, size_t fromCol) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), 0, toRow, fromCol, self.getRow() - fromCol);
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::bottomLeftCorner(this auto&& self, size_t fromRow, size_t toCol) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), fromRow, self.getRow() - fromRow, 0, toCol);
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::bottomRightCorner(this auto&& self, size_t fromRow, size_t fromCol) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), fromRow, self.getRow() - fromRow, fromCol, self.getCol() - fromCol);
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::bottomRightCorner(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), from, self.getRow() - from, from, self.getCol() - from);
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::block(this auto&& self, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::diag(this auto&& self) noexcept {
        using Self = decltype(self);
        return MainDiag<Self>(std::forward<Self>(self));
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::diag(this auto&& self, ssize_t offset) noexcept {
        using Self = decltype(self);
        return OffsetDiag<Self>(std::forward<Self>(self), offset);
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::triu(this auto&& self) noexcept {
        using Self = decltype(self);
        return MatrixTrig<Self, true, false>(std::forward<Self>(self));
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::triu_unit(this auto&& self) noexcept {
        using Self = decltype(self);
        return MatrixTrig<Self, true, true>(std::forward<Self>(self));
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::tril(this auto&& self) noexcept {
        using Self = decltype(self);
        return MatrixTrig<Self, false, false>(std::forward<Self>(self));
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::tril_unit(this auto&& self) noexcept {
        using Self = decltype(self);
        return MatrixTrig<Self, false, true>(std::forward<Self>(self));
    }

    template<class Derived, Scalar ScalarT>
    Index2D RValueMatrix<Derived, ScalarT>::argmax() const noexcept {
        Trv x = std::numeric_limits<Trv>::lowest();
        Index2D result{0, 0};
        for (size_t r = 0; r < getRow(); ++r) {
            for (size_t c = 0; c < getCol(); ++c) {
                Trv y = calc_value(r, c);
                if (y > x) {
                    x = y;
                    result = Index2D{r, c};
                }
            }
        }
        return result;
    }

    template<class Derived, Scalar ScalarT>
    Index2D RValueMatrix<Derived, ScalarT>::argmin() const noexcept {
        Trv x = std::numeric_limits<Trv>::max();
        Index2D result{0, 0};
        for (size_t r = 0; r < getRow(); ++r) {
            for (size_t c = 0; c < getCol(); ++c) {
                Trv y = calc_value(r, c);
                if (y < x) {
                    x = y;
                    result = Index2D{r, c};
                }
            }
        }
        return result;
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::max() const -> T {
        T result;
        if constexpr (Derived::isColMatrix()) {
            result = Base::getDerived().col(0).max();
            for (size_t i = 1; i < getCol(); ++i) {
                T temp = Base::getDerived().col(i).max();
                if (temp > result)
                    result = temp;
            }
        }
        else {
            result = Base::getDerived().row(0).max();
            for (size_t i = 1; i < getRow(); ++i) {
                T temp = Base::getDerived().row(i).max();
                if (temp > result)
                    result = temp;
            }
        }
        return result;
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::min() const -> T {
        T result;
        if constexpr (Derived::isColMatrix()) {
            result = Base::getDerived().col(0).min();
            for (size_t i = 1; i < getCol(); ++i) {
                T temp = Base::getDerived().col(i).min();
                if (temp < result)
                    result = temp;
            }
        }
        else {
            result = row(0).min();
            for (size_t i = 1; i < getRow(); ++i) {
                T temp = Base::getDerived().row(i).min();
                if (temp < result)
                    result = temp;
            }
        }
        return result;
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::sum() const -> CoDiff<T> {
        const auto& x = Base::getDerived();
        if constexpr (isReverseDiff()) {
            auto& result = co_yield x.values().sum();
            x.reverse(result.grad());
        }
        else {
            T result = 0;
            for (size_t major = 0; major < getMaxMajor(); ++major)
                result += MatrixMajor::isColMatrix<Derived>() ? x.col(major).sum() : x.row(major).sum();
            co_return std::move(result);
        }
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::sum_rows() const {
        return MatrixSum<Derived, false>(Base::getDerived());
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::sum_cols() const {
        return MatrixSum<Derived, true>(Base::getDerived());
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::mean() const -> CoDiff<T> {
        return sum() / Trv(getSize());
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::trace() const -> T {
        assert(isSquare());
        T result = T(0);
        for (size_t i = 0; i < getRow(); ++i)
            result += calc(i, i);
        return result;
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::lnSumExp() const -> CoDiff<T> {
        return Base::getDerived().flatten().lnSumExp();
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::det() const -> CoDiff<T> {
        assert(isSquare() && "[Error]: Determinate requires square matrix");
        constexpr size_t Order = std::max(Derived::getRowAtCompile(), Derived::getColAtCompile());
        if constexpr (Order == 1)
            return calc(0, 0);
        else if constexpr (Order == 2)
            return fma(calc(0, 0), calc(1, 1), -calc(0, 1) * calc(1, 0));
        else if constexpr (Order == 3)
            return fma(calc(0, 0), fma(calc(1, 1), calc(2, 2), -calc(1, 2) * calc(2, 1)), fma(calc(0, 1), fma(calc(1, 2), calc(2, 0), -calc(1, 0) * calc(2, 2)), calc(0, 2) * fma(calc(1, 0), calc(2, 1), -calc(1, 1) * calc(2, 0))));
        else
            return DenseLU<T, false>(Base::getDerived()).det();
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::lnAbsDet() const -> Tr {
        constexpr size_t Order = std::max(Derived::getRowAtCompile(), Derived::getColAtCompile());
        if constexpr (Order <= 3)
            return ln(abs(Base::getDerived().det()));
        else {
            DenseLU<T, false> lu(Base::getDerived());
            return ln(abs(lu.getMatrixLU().diag())).sum();
        }
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::sgndet() const {
        DenseLU<T, false> lu(Base::getDerived());
        return unit(lu.getMatrixLU().diag()).prod();
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::format() const noexcept {
        return FormatedMatrix<Derived>(Base::getDerived());
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::inv(this auto&& self) noexcept {
        using Self = decltype(self);
        return Inverse<Self>(std::forward<Self>(self));
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::pinv() const noexcept {
        return PseudoInverse<Derived>(Base::getDerived());
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::transpose(this auto&& self) noexcept {
        using Self = decltype(self);
        return Transpose<Self>(std::forward<Self>(self));
    }

    template<class Derived, Scalar ScalarT>
    decltype(auto) RValueMatrix<Derived, ScalarT>::conjugate(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (isComplex())
            return Conjugate<Self>(std::forward<Self>(self));
        else
            return std::forward<Self>(self);
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::hermite(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (isComplex())
            return Hermite<Self>(std::forward<Self>(self));
        else
            return std::forward<Self>(self).transpose();
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::flatten(this auto&& self) noexcept {
        using Self = decltype(self);
        return Flatten<Self>(std::forward<Self>(self));
    }

    template<class Derived, Scalar ScalarT>
    decltype(auto) RValueMatrix<Derived, ScalarT>::reals(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (isComplex())
            return RealMatrix<Self>(std::forward<Self>(self));
        else
            return std::forward<Self>(self);
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::imags(this auto&& self) noexcept {
        using Self = decltype(self);
        return ImagMatrix<Self>(std::forward<Self>(self));
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::squaredNorms(this auto&& self) noexcept {
        using Self = decltype(self);
        return SquaredNormMatrix<Self>(std::forward<Self>(self));
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::norms(this auto&& self) noexcept {
        using Self = decltype(self);
        return NormMatrix<Self>(std::forward<Self>(self));
    }

    template<class Derived, Scalar ScalarT>
    decltype(auto) RValueMatrix<Derived, ScalarT>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (isDiffable())
            return ValueMatrix<Self>(std::forward<Self>(self));
        else
            return std::forward<Self>(self);
    }

    template<class Derived, Scalar ScalarT>
    template<int GradOrder>
    decltype(auto) RValueMatrix<Derived, ScalarT>::grads(this auto&& self) noexcept {
        using Self = decltype(self);
        return GradMatrix<Self, GradOrder>(std::forward<Self>(self));
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::getRow() const noexcept {
        if constexpr (Derived::isStaticSquare())
            return getOrder();
        else
            return Base::getDerived().getRow();
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::getCol() const noexcept {
        if constexpr (Derived::isStaticSquare())
            return getOrder();
        else
            return Base::getDerived().getCol();
    }

    template<class Derived, Scalar ScalarT>
    auto RValueMatrix<Derived, ScalarT>::getOrder() const noexcept {
        assert(isSquare() && "[Error]: getOrder() assumes square matrix");
        return Base::getDerived().getOrder();
    }

    template<class Derived, Scalar ScalarT>
    size_t RValueMatrix<Derived, ScalarT>::getMaxMajor(this const auto& self) noexcept {
        return Derived::isColMatrix() ? self.getCol() : self.getRow();
    }

    template<class Derived, Scalar ScalarT>
    size_t RValueMatrix<Derived, ScalarT>::getMaxMinor(this const auto& self) noexcept {
        return Derived::isColMatrix() ? self.getRow() : self.getCol();
    }

    template<class Derived, Scalar ScalarT>
    constexpr bool RValueMatrix<Derived, ScalarT>::isOverdetermined() const noexcept {
        if constexpr (Derived::isStaticSquare())
            return false;
        else
            return getRow() > getCol();
    }

    template<class Derived, Scalar ScalarT>
    constexpr bool RValueMatrix<Derived, ScalarT>::isUnderdetermined() const noexcept {
        if constexpr (Derived::isStaticSquare())
            return false;
        else
            return getRow() < getCol();
    }

    template<class Derived, Scalar ScalarT>
    constexpr bool RValueMatrix<Derived, ScalarT>::isSquare() const noexcept {
        if constexpr (Derived::isStaticSquare())
            return true;
        else
            return getRow() == getCol();
    }

    template<class Derived, Scalar ScalarT>
    constexpr bool RValueMatrix<Derived, ScalarT>::isSymm() const noexcept {
        if constexpr (Derived::isStaticSymm())
            return true;
        else {
            if (!isSquare())
                return false;

            const size_t order = getRow();
            for (size_t r = 0; r < order; ++r)
                for (size_t c = r + 1; c < order; ++c)
                    if (calc(r, c) != calc(c, r))
                        return false;
            return true;
        }
    }

    template<class Derived, Scalar ScalarT>
    constexpr bool RValueMatrix<Derived, ScalarT>::isHermite() const noexcept {
        if constexpr (Derived::isStaticHermite())
            return true;
        else {
            if (!isSquare())
                return false;

            const size_t order = getRow();
            for (size_t r = 0; r < order; ++r) {
                if (!calc(r, r).imag().isZero())
                    return false;

                for (size_t c = r + 1; c < order; ++c)
                    if (calc(r, c) != calc(c, r).conjugate())
                        return false;
            }
            return true;
        }
    }

    template<class Derived, Scalar ScalarT>
    bool RValueMatrix<Derived, ScalarT>::isFinite() const noexcept {
        for (size_t r = 0; r < getRow(); ++r)
            for (size_t c = 0; c < getCol(); ++c)
                if (!calc(r, c).isFinite())
                    return false;
        return true;
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueMatrix<Derived, ScalarT>::isComplex() noexcept {
        return ScalarType::isComplex();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueMatrix<Derived, ScalarT>::isDiffable() noexcept {
        return ScalarType::isDiffable();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueMatrix<Derived, ScalarT>::isForwardDiff() noexcept {
        return ScalarType::isForwardDiff();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueMatrix<Derived, ScalarT>::isReverseDiff() noexcept {
        return ScalarType::isReverseDiff();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueMatrix<Derived, ScalarT>::isLValueMatrix() noexcept {
        return false;
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueMatrix<Derived, ScalarT>::isStrided() noexcept {
        return false;
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueMatrix<Derived, ScalarT>::isCompact() noexcept {
        return false;
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueMatrix<Derived, ScalarT>::isSparse() noexcept {
        return requires { std::declval<Derived>().getNumNonzero(); };
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueMatrix<Derived, ScalarT>::isStaticSymm() noexcept {
        using TransposeTy = std::remove_cvref_t<decltype(std::declval<const Derived&>().transpose())>;
        return std::same_as<TransposeTy, Derived>;
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueMatrix<Derived, ScalarT>::isStaticHermite() noexcept {
        if constexpr (T::isComplex()) {
            using HermiteTy = std::remove_cvref_t<decltype(std::declval<const Derived&>().hermite())>;
            return std::same_as<HermiteTy, Derived>;
        }
        else
            return isStaticSymm();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueMatrix<Derived, ScalarT>::isStaticSquare() noexcept {
        if constexpr (isStaticHermite())
            return true;
        else
            return Derived::getRowAtCompile() != Dynamic && Derived::getRowAtCompile() == Derived::getColAtCompile();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueMatrix<Derived, ScalarT>::isColMatrix() noexcept {
        return (Derived::getMajor() & MatrixMajor::Col) != 0;
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueMatrix<Derived, ScalarT>::isRowMatrix() noexcept {
        return (Derived::getMajor() & MatrixMajor::Row) != 0;
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueMatrix<Derived, ScalarT>::isBothMajor() noexcept {
        return Derived::getMajor() == MatrixMajor::BothMajor;
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool RValueMatrix<Derived, ScalarT>::isFastAssign() noexcept {
        return false;
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval size_t RValueMatrix<Derived, ScalarT>::getRowAtCompile() noexcept {
        return Dynamic;
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval size_t RValueMatrix<Derived, ScalarT>::getColAtCompile() noexcept {
        return Dynamic;
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval size_t RValueMatrix<Derived, ScalarT>::getSizeAtCompile() noexcept {
        return Derived::getRowAtCompile() * Derived::getColAtCompile();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval auto RValueMatrix<Derived, ScalarT>::getMajor() noexcept {
        return Derived::getMajor();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval void RValueMatrix<Derived, ScalarT>::static_assert_assign(const Scalar auto& source) noexcept {
        using U = std::remove_cvref<decltype(source)>::type;
        T::template static_assert_assign<U>();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval void RValueMatrix<Derived, ScalarT>::static_assert_assign(const Matrix auto& source) noexcept {
        constexpr size_t R1 = Derived::getRowAtCompile();
        constexpr size_t C1 = Derived::getColAtCompile();
        constexpr size_t R2 = source.getRowAtCompile();
        constexpr size_t C2 = source.getColAtCompile();
        static_assert(R1 == R2 || R1 == Dynamic || R2 == Dynamic, "[Error]: Row mismatch between two matrix");
        static_assert(C1 == C2 || C1 == Dynamic || C2 == Dynamic, "[Error]: Col mismatch between two matrix");

        using U = std::remove_cvref_t<decltype(source)>::ScalarType;
        T::template static_assert_assign<U>();
    }

    template<class Derived, Scalar ScalarT>
    consteval int RValueMatrix<Derived, ScalarT>::calcBlockingSize(int CacheSize) noexcept {
        int result = 1;
        while (result * result * int(sizeof(Trv)) < CacheSize)
            result *= 2;
        return result / 2;
    }
    /**
     * See if the block range is legal to the matrix
     */
    template<class Derived, Scalar ScalarT>
    __host__ __device__ void RValueMatrix<Derived, ScalarT>::checkBlock(
            [[maybe_unused]] const Matrix auto& m,
            [[maybe_unused]] size_t fromRow,
            [[maybe_unused]] size_t rowCount,
            [[maybe_unused]] size_t fromCol,
            [[maybe_unused]] size_t colCount) noexcept {
        assert(fromRow < m.getRow());
        assert(fromCol < m.getCol());
        assert((fromRow + rowCount) <= m.getRow());
        assert((fromCol + colCount) <= m.getCol());
    }
}

#include "Sum.h"
#include "Inverse.h"
#include "PseudoInverse.h"
#include "Transpose.h"
#include "Conjugate.h"
#include "Hermite.h"
#include "Flatten.h"
#include "Trig/MatrixTrig.h"
#include "Convert/MatrixConvert.h"
#include "ReshapedVector.h"
