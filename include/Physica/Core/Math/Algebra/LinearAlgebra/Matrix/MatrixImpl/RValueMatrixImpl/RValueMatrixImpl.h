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
    template<class Derived>
    auto RValueMatrix<Derived>::operator*(this auto&& self, Vector auto&& v) noexcept {
        assert(self.getCol() == v.getLength());
        using Self = decltype(self);
        using V = decltype(v);
        static_assert(!is_device_obj<V>::value, "[Error]: host-device mismatch");
        if constexpr (RowAtCompile == 1)
            return std::forward<Self>(self).row(0) * v;
        else
            return GEMV<Self, V&&>(std::forward<Self>(self), std::forward<V>(v));
    }

    template<class Derived>
    auto RValueMatrix<Derived>::operator*(this auto&& self, Matrix auto&& m) noexcept {
        using Self = decltype(self);
        using M = decltype(m);
        static_assert(!is_device_obj<M>::value, "[Error]: host-device mismatch");
        constexpr bool ColVectorLHS = ColAtCompile == 1;
        constexpr bool RowVectorLHS = RowAtCompile == 1;
        constexpr bool ColVectorRHS = Traits<M>::ColAtCompile == 1;
        constexpr bool RowVectorRHS = Traits<M>::RowAtCompile == 1;
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

    template<class Derived>
    template<ExecutePolicy P>
    void RValueMatrix<Derived>::assign(Matrix auto&& target) const noexcept {
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

    template<class Derived>
    template<ExecutePolicy P>
    void RValueMatrix<Derived>::assign_base(Matrix auto&& target) const noexcept {
        assign<P>(target);
    }

    template<class Derived>
    template<ExecutePolicy P>
    void RValueMatrix<Derived>::assign_add(Matrix auto&& target) const noexcept {
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

    template<class Derived>
    void RValueMatrix<Derived>::assert_assign(const Matrix auto& source) const noexcept {
        static_assert_assign(source);

        constexpr size_t Row1 = RowAtCompile;
        constexpr size_t Row2 = source.RowAtCompile;
        if constexpr (Row1 == Dynamic || Row2 == Dynamic)
            assert(getRow() == source.getRow() && "[Error]: Dimensions do not match");

        constexpr size_t Col1 = RowAtCompile;
        constexpr size_t Col2 = source.RowAtCompile;
        if constexpr (Col1 == Dynamic || Col2 == Dynamic)
            assert(getCol() == source.getCol() && "[Error]: Dimensions do not match");
        
        constexpr size_t Size1 = SizeAtCompile;
        constexpr size_t Size2 = source.SizeAtCompile;
        if constexpr (Size1 == Dynamic || Size2 == Dynamic)
            assert(getSize() > 0);
    }

    template<class Derived>
    void RValueMatrix<Derived>::assert_assign_mkl(const Matrix auto& source) const noexcept {
        static_assert(Internal::EnableMKL<Derived, decltype(source)>::value, "[Error]: Cannot apply MKL to this expr");
        assert_assign(source);
    }

    template<class Derived>
    decltype(auto) RValueMatrix<Derived>::calc(size_t row, size_t col) const {
        return Base::getDerived().calc(row, col);
    }

    template<class Derived>
    decltype(auto) RValueMatrix<Derived>::calc_value(size_t row, size_t col) const {
        return Base::getDerived().values().calc(row, col);
    }

    template<class Derived>
    decltype(auto) RValueMatrix<Derived>::calcFromMajorMinor(size_t major, size_t minor) const {
        return calc(rowFromMajorMinor(major, minor), colFromMajorMinor(major, minor));
    }

    template<class Derived>
    void RValueMatrix<Derived>::reverse(const Matrix auto&, const Matrix auto& grad) const noexcept {
        static_assert(isReverseDiff());
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
    auto RValueMatrix<Derived>::row(this auto&& self, size_t r) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self, 1, Dynamic>(std::forward<Self>(self), r, 0, self.getCol());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::col(this auto&& self, size_t c) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self, Dynamic, 1>(std::forward<Self>(self), 0, self.getRow(), c);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::rows(this auto&& self, size_t fromRow, size_t rowCount) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), fromRow, rowCount, 0, self.getCol());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::topRows(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), 0, to, 0, self.getCol());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::bottomRows(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), from, self.getRow() - from, 0, self.getCol());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::cols(this auto&& self, size_t fromCol, size_t colCount) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), 0, self.getRow(), fromCol, colCount);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::leftCols(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), 0, self.getRow(), 0, to);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::rightCols(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), 0, self.getRow(), from, self.getCol() - from);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::topLeftCorner(this auto&& self, size_t toRow, size_t toCol) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), 0, toRow, 0, toCol);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::topLeftCorner(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), 0, to, 0, to);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::topRightCorner(this auto&& self, size_t toRow, size_t fromCol) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), 0, toRow, fromCol, self.getRow() - fromCol);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::bottomLeftCorner(this auto&& self, size_t fromRow, size_t toCol) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), fromRow, self.getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::bottomRightCorner(this auto&& self, size_t fromRow, size_t fromCol) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), fromRow, self.getRow() - fromRow, fromCol, self.getCol() - fromCol);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::bottomRightCorner(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), from, self.getRow() - from, from, self.getCol() - from);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::block(this auto&& self, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        using Self = decltype(self);
        return RMatrixBlock<Self>(std::forward<Self>(self), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::diag(this auto&& self) noexcept {
        assert(self.isSquare());
        using Self = decltype(self);
        return DiagVectorR<Self>(std::forward<Self>(self));
    }

    template<class Derived>
    auto RValueMatrix<Derived>::diag(this auto&& self, ssize_t shift) noexcept {
        using Self = decltype(self);
        return MinorDiagR<Self>(std::forward<Self>(self), shift);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::triu(this auto&& self) noexcept {
        using Self = decltype(self);
        return MatrixTrig<Self, true, false>(std::forward<Self>(self));
    }

    template<class Derived>
    auto RValueMatrix<Derived>::triu_unit(this auto&& self) noexcept {
        using Self = decltype(self);
        return MatrixTrig<Self, true, true>(std::forward<Self>(self));
    }

    template<class Derived>
    auto RValueMatrix<Derived>::tril(this auto&& self) noexcept {
        using Self = decltype(self);
        return MatrixTrig<Self, false, false>(std::forward<Self>(self));
    }

    template<class Derived>
    auto RValueMatrix<Derived>::tril_unit(this auto&& self) noexcept {
        using Self = decltype(self);
        return MatrixTrig<Self, false, true>(std::forward<Self>(self));
    }

    template<class Derived>
    Index2D RValueMatrix<Derived>::argmax() const noexcept {
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

    template<class Derived>
    Index2D RValueMatrix<Derived>::argmin() const noexcept {
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

    template<class Derived>
    auto RValueMatrix<Derived>::max() const -> T {
        T result;
        if constexpr (MatrixMajor::isColMatrix<This>()) {
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

    template<class Derived>
    auto RValueMatrix<Derived>::min() const -> T {
        T result;
        if constexpr (MatrixMajor::isColMatrix<This>()) {
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

    template<class Derived>
    auto RValueMatrix<Derived>::sum() const -> CoDiff<T> {
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
        return Base::getDerived().flatten().lnSumExp();
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
            return DenseLU<T, false>(Base::getDerived()).det();
    }

    template<class Derived>
    auto RValueMatrix<Derived>::lnAbsDet() const -> Tr {
        DenseLU<T, false> lu(Base::getDerived());
        return ln(abs(lu.getMatrixLU().diag())).sum();
    }

    template<class Derived>
    auto RValueMatrix<Derived>::sgndet() const {
        DenseLU<T, false> lu(Base::getDerived());
        return unit(lu.getMatrixLU().diag()).prod();
    }

    template<class Derived>
    auto RValueMatrix<Derived>::format() const noexcept {
        return FormatedMatrix<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::inv(this auto&& self) noexcept {
        using Self = decltype(self);
        return Inverse<Self>(std::forward<Self>(self));
    }

    template<class Derived>
    auto RValueMatrix<Derived>::pinv() const noexcept {
        return PseudoInverse<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueMatrix<Derived>::transpose(this auto&& self) noexcept {
        using Self = decltype(self);
        return Transpose<Self>(std::forward<Self>(self));
    }

    template<class Derived>
    decltype(auto) RValueMatrix<Derived>::conjugate(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (isComplex())
            return Conjugate<Self>(std::forward<Self>(self));
        else
            return std::forward<Self>(self);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::hermite(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (isComplex())
            return Hermite<Self>(std::forward<Self>(self));
        else
            return std::forward<Self>(self).transpose();
    }

    template<class Derived>
    auto RValueMatrix<Derived>::flatten(this auto&& self) noexcept {
        using Self = decltype(self);
        return FlattenR<Self>(std::forward<Self>(self));
    }

    template<class Derived>
    decltype(auto) RValueMatrix<Derived>::reals(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (isComplex())
            return RealMatrix<Self>(std::forward<Self>(self));
        else
            return std::forward<Self>(self);
    }

    template<class Derived>
    auto RValueMatrix<Derived>::imags(this auto&& self) noexcept {
        using Self = decltype(self);
        return ImagMatrix<Self>(std::forward<Self>(self));
    }
    
    template<class Derived>
    auto RValueMatrix<Derived>::squaredNorms(this auto&& self) noexcept {
        using Self = decltype(self);
        return SquaredNormMatrix<Self>(std::forward<Self>(self));
    }

    template<class Derived>
    auto RValueMatrix<Derived>::norms(this auto&& self) noexcept {
        using Self = decltype(self);
        return NormMatrix<Self>(std::forward<Self>(self));
    }

    template<class Derived>
    decltype(auto) RValueMatrix<Derived>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (isDiffable())
            return ValueMatrix<Self>(std::forward<Self>(self));
        else
            return std::forward<Self>(self);
    }

    template<class Derived>
    template<int GradOrder>
    decltype(auto) RValueMatrix<Derived>::grads(this auto&& self) noexcept {
        using Self = decltype(self);
        return GradMatrix<Self, GradOrder>(std::forward<Self>(self));
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
        if constexpr (isStaticSymm())
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

    template<class Derived>
    bool RValueMatrix<Derived>::isHermite() const noexcept {
        if constexpr (isStaticHermite())
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

    template<class Derived>
    bool RValueMatrix<Derived>::isFinite() const noexcept {
        for (size_t r = 0; r < getRow(); ++r)
            for (size_t c = 0; c < getCol(); ++c)
                if (!calc(r, c).isFinite())
                    return false;
        return true;
    }

    template<class Derived>
    __host__ __device__ consteval bool RValueMatrix<Derived>::isComplex() noexcept {
        return ScalarType::isComplex();
    }

    template<class Derived>
    __host__ __device__ consteval bool RValueMatrix<Derived>::isDiffable() noexcept {
        return ScalarType::isDiffable();
    }

    template<class Derived>
    __host__ __device__ consteval bool RValueMatrix<Derived>::isForwardDiff() noexcept {
        return ScalarType::isForwardDiff();
    }

    template<class Derived>
    __host__ __device__ consteval bool RValueMatrix<Derived>::isReverseDiff() noexcept {
        return ScalarType::isReverseDiff();
    }

    template<class Derived>
    __host__ __device__ consteval bool RValueMatrix<Derived>::isLValueMatrix() noexcept {
        return false;
    }

    template<class Derived>
    __host__ __device__ consteval bool RValueMatrix<Derived>::isCompact() noexcept {
        return requires{ std::declval<Derived>().data(); };
    }

    template<class Derived>
    __host__ __device__ consteval bool RValueMatrix<Derived>::isSparse() noexcept {
        return requires{ std::declval<Derived>().getNumNonzero(); };
    }

    template<class Derived>
    __host__ __device__ consteval bool RValueMatrix<Derived>::isStaticSymm() noexcept {
        using TransposeTy = std::remove_cvref_t<decltype(std::declval<Derived>().transpose())>;
        return std::same_as<TransposeTy, Derived>;
    }

    template<class Derived>
    __host__ __device__ consteval bool RValueMatrix<Derived>::isStaticHermite() noexcept {
        using HermiteTy = std::remove_cvref_t<decltype(std::declval<Derived>().hermite())>;
        return std::same_as<HermiteTy, Derived>;
    }

    template<class Derived>
    __host__ __device__ consteval bool RValueMatrix<Derived>::isColMatrix() noexcept {
        return (getMajor() & MatrixMajor::Col) != 0;
    }

    template<class Derived>
    __host__ __device__ consteval bool RValueMatrix<Derived>::isRowMatrix() noexcept {
        return (getMajor() & MatrixMajor::Row) != 0;
    }

    template<class Derived>
    __host__ __device__ consteval bool RValueMatrix<Derived>::isBothMajor() noexcept {
        return getMajor() == MatrixMajor::BothMajor;
    }

    template<class Derived>
    __host__ __device__ consteval int RValueMatrix<Derived>::getMajor() noexcept {
        return Traits<Derived>::Major;
    }

    template<class Derived>
    __host__ __device__ consteval void RValueMatrix<Derived>::static_assert_assign(const Scalar auto& source) noexcept {
        using U = std::remove_cvref<decltype(source)>::type;
        T::template static_assert_assign<U>();
    }

    template<class Derived>
    __host__ __device__ consteval void RValueMatrix<Derived>::static_assert_assign(const Matrix auto& source) noexcept {
        using Src = std::remove_cvref_t<decltype(source)>;
        static_assert(RowAtCompile == Src::RowAtCompile || RowAtCompile == Dynamic || Src::RowAtCompile == Dynamic, "[Error]: Row mismatch between two matrix");
        static_assert(ColAtCompile == Src::ColAtCompile || ColAtCompile == Dynamic || Src::ColAtCompile == Dynamic, "[Error]: Col mismatch between two matrix");

        using U = Src::ScalarType;
        T::template static_assert_assign<U>();
    }

    template<class Derived>
    consteval int RValueMatrix<Derived>::calcBlockingSize(int CacheSize) noexcept {
        int result = 1;
        while (result * result * int(sizeof(Trv)) < CacheSize)
            result *= 2;
        return result / 2;
    }
    /**
     * See if the block range is legal to the matrix
     */
     template<class Derived>
    __host__ __device__ void RValueMatrix<Derived>::checkBlock(
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
