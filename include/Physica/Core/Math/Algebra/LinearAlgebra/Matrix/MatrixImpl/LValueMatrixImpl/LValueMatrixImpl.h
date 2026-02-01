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

#include "../LValueMatrix.h"
#include "Flatten.h"

namespace Physica {
    template<class Derived>
    auto LValueMatrix<Derived>::operator=(Scalar auto x) noexcept -> Derived& {
        Base::static_assert_assign(x);
        const size_t maxMajor = Base::getMaxMajor();
        const size_t maxMinor = Base::getMaxMinor();
        for (size_t i = 0; i < maxMajor; ++i)
            for (size_t j = 0; j < maxMinor; ++j)
                    refFromMajorMinor(i, j) = x;
        return Base::getDerived();
    }

    template<class Derived>
    void LValueMatrix<Derived>::operator+=(Scalar auto x) noexcept {
        auto& m = Base::getDerived();
        (m + x).assign(m);
    }

    template<class Derived>
    void LValueMatrix<Derived>::operator-=(Scalar auto x) noexcept {
        auto& m = Base::getDerived();
        (m - x).assign(m);
    }

    template<class Derived>
    void LValueMatrix<Derived>::operator*=(Scalar auto x) noexcept {
        auto& m = Base::getDerived();
        (m * x).assign(m);
    }

    template<class Derived>
    void LValueMatrix<Derived>::operator/=(Scalar auto x) noexcept {
        auto& m = Base::getDerived();
        (m / x).assign(m);
    }

    template<class Derived>
    Derived& LValueMatrix<Derived>::operator=(const Matrix auto& m) {
        if constexpr (std::is_same<const Derived&, decltype(m)>::value)
            assert(this != &m && "[Error]: Self assign is likely a bug");
        Base::getDerived().resize(m);
        m.assign(Base::getDerived());
        return Base::getDerived();
    }

    template<class Derived>
    void LValueMatrix<Derived>::operator+=(const Matrix auto& m) {
        m.assign_add(Base::getDerived());
    }

    template<class Derived>
    void LValueMatrix<Derived>::operator-=(const Matrix auto& m) {
        Base::getDerived() += -m;
    }

    template<class Derived>
    decltype(auto) LValueMatrix<Derived>::operator[](this auto&& self, size_t row, size_t col) {
        return *self.data_ptr(row, col);
    }

    template<class Derived>
    size_t LValueMatrix<Derived>::pivotPartial(size_t column) const noexcept {
        assert(column < Base::getCol());
        assert(Base::isSquare());
        return column + abs(Base::getDerived().col(column).tail(column)).argmax();
    }

    template<class Derived>
    size_t LValueMatrix<Derived>::pivotComplete(size_t column) const noexcept {
        assert(column < Base::getCol());
        assert(Base::isSquare());
        return column + abs_elem(Base::getDerived().bottomRightCorner(column)).argmax()[0];
    }

    template<class Derived>
    auto LValueMatrix<Derived>::sum() const -> CoDiff<ScalarType> {
        if constexpr (isReverseDiff) {
            const size_t maxMajor = Base::getMaxMajor();
            const size_t maxMinor = Base::getMaxMinor();
            Tr v = 0;
            for (size_t major = 0; major < maxMajor; ++major)
                for (size_t minor = 0; minor < maxMinor; ++minor)
                    v += refFromMajorMinor(major, minor).value();

            const auto result = co_yield std::move(v);
            for (size_t major = 0; major < maxMajor; ++major)
                for (size_t minor = 0; minor < maxMinor; ++minor)
                    refFromMajorMinor(major, minor).reverse(result.grad());
        }
        else
            co_return Base::sum();
    }

    template<class Derived>
    void LValueMatrix<Derived>::toNextMean(size_t lastNumSample, const Matrix auto& sample) noexcept {
        Base::assert_assign(sample);
        auto& mean = Base::getDerived();
        if (lastNumSample == 0) [[unlikely]] {
            sample.assign(mean);
            return;
        }

        const T factor1 = T(lastNumSample);
        const T factor2 = reciprocal(T(lastNumSample + 1));
        mean = (factor1 * mean + sample) * factor2;
    }

    template<class Derived>
    void LValueMatrix<Derived>::toNextVariance(Derived& mean, size_t lastNumSample, const Matrix auto& sample) noexcept {
        Base::assert_assign(sample);
        auto& var = Base::getDerived();
        if (lastNumSample == 0) [[unlikely]] {
            sample.assign(mean);
            var.zeros();
            return;
        }

        const T factor1 = T(lastNumSample);
        const T factor2 = reciprocal(T(lastNumSample + 1));
        var = (var + square_elem(mean - sample) * factor2) * (factor1 * factor2);
        mean.toNextMean(lastNumSample, sample);
    }

    template<class Derived>
    void LValueMatrix<Derived>::reverse(const auto& grad) const noexcept {
        using U = std::remove_cvref_t<decltype(grad)>;
        static_assert(std::same_as<typename ScalarType::GradType, typename U::ScalarType>, "[Error]: Inconsistent ScalarType");
        static_assert(isReverseDiff);
        if constexpr (Scalar<U>)
            Base::getConstCastDerived().grads() += grad;
        else {
            static_assert(Matrix<U>, "[Error]: Unexpected type");
            Base::assert_assign(grad);
            grad.assign_add(Base::getConstCastDerived().grads());
        }
    }

    template<class Derived>
    auto LValueMatrix<Derived>::row(this auto&& self, size_t r) noexcept {
        using Self = decltype(self);
        return LMatrixBlock<Self, 1, Dynamic>(std::forward<Self>(self), r, 0, self.getCol());
    }

    template<class Derived>
    auto LValueMatrix<Derived>::col(this auto&& self, size_t c) noexcept {
        using Self = decltype(self);
        return LMatrixBlock<Self, Dynamic, 1>(std::forward<Self>(self), 0, self.getRow(), c);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::rows(this auto&& self, size_t fromRow, size_t rowCount) noexcept {
        using Self = decltype(self);
        return LMatrixBlock<Self>(std::forward<Self>(self), fromRow, rowCount, 0, self.getCol());
    }

    template<class Derived>
    auto LValueMatrix<Derived>::topRows(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return LMatrixBlock<Self>(std::forward<Self>(self), 0, to, 0, self.getCol());
    }

    template<class Derived>
    auto LValueMatrix<Derived>::bottomRows(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return LMatrixBlock<Self>(std::forward<Self>(self), from, self.getRow() - from, 0, self.getCol());
    }

    template<class Derived>
    auto LValueMatrix<Derived>::cols(this auto&& self, size_t fromCol, size_t colCount) noexcept {
        using Self = decltype(self);
        return LMatrixBlock<Self>(std::forward<Self>(self), 0, self.getRow(), fromCol, colCount);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::leftCols(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return LMatrixBlock<Self>(std::forward<Self>(self), 0, self.getRow(), 0, to);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::rightCols(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return LMatrixBlock<Self>(std::forward<Self>(self), 0, self.getRow(), from, self.getCol() - from);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::topLeftCorner(this auto&& self, size_t toRow, size_t toCol) noexcept {
        using Self = decltype(self);
        return LMatrixBlock<Self>(std::forward<Self>(self), 0, toRow, 0, toCol);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::topLeftCorner(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return LMatrixBlock<Self>(std::forward<Self>(self), 0, to, 0, to);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::topRightCorner(this auto&& self, size_t toRow, size_t fromCol) noexcept {
        using Self = decltype(self);
        return LMatrixBlock<Self>(std::forward<Self>(self), 0, toRow, fromCol, self.getRow() - fromCol);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::bottomLeftCorner(this auto&& self, size_t fromRow, size_t toCol) noexcept {
        using Self = decltype(self);
        return LMatrixBlock<Self>(std::forward<Self>(self), fromRow, self.getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::bottomRightCorner(this auto&& self, size_t fromRow, size_t fromCol) noexcept {
        using Self = decltype(self);
        return LMatrixBlock<Self>(std::forward<Self>(self), fromRow, self.getRow() - fromRow, fromCol, self.getCol() - fromCol);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::bottomRightCorner(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return LMatrixBlock<Self>(std::forward<Self>(self), from, self.getRow() - from, from, self.getCol() - from);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::block(this auto&& self, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        using Self = decltype(self);
        return LMatrixBlock<Self>(std::forward<Self>(self), fromRow, rowCount, fromCol, colCount);
    }


    template<class Derived>
    auto LValueMatrix<Derived>::diag(this auto&& self) noexcept {
        assert(self.isSquare());
        using Self = decltype(self);
        return DiagVectorL<Self>(std::forward<Self>(self));
    }

    template<class Derived>
    auto LValueMatrix<Derived>::diag(this auto&& self, ssize_t shift) noexcept {
        using Self = decltype(self);
        return MinorDiagL<Self>(std::forward<Self>(self), shift);
    }
    /**
     * Reduce the element at one row using the other row.
     * \param r1
     * The index of row to be used.
     * \param r2
     * The index of row that the element belongs to.
     * \param elementIndex
     * Index of the element to be reduced.
     */
    template<class Derived>
    void LValueMatrix<Derived>::rowReduce(size_t r1, size_t r2, size_t elementIndex) {
        Derived& matrix = Base::getDerived();
        assert(abs(matrix[r1, elementIndex]) > ScalarType(std::numeric_limits<ScalarType>::min()));
        const ScalarType factor = matrix[r2, elementIndex] / matrix[r1, elementIndex];
        matrix.row(r2) -= matrix.row(r1) * factor;
        matrix[r2, elementIndex] = 0;
    }

    template<class Derived>
    void LValueMatrix<Derived>::colReduce(size_t c1, size_t c2, size_t elementIndex) {
        Derived& matrix = Base::getDerived();
        assert(abs(matrix(elementIndex, c1)) > ScalarType(std::numeric_limits<ScalarType>::min()));
        const ScalarType factor = matrix(c2, elementIndex) / matrix(c1, elementIndex);
        matrix.col(c2) -= matrix.col(c1) * factor;
        matrix[elementIndex, c2] = 0;
    }

    template<class Derived>
    void LValueMatrix<Derived>::majorReduce(size_t v1, size_t v2, size_t elementIndex) {
        if constexpr (MatrixOption::isColMatrix<Derived>())
            colReduce(v1, v2, elementIndex);
        else
            rowReduce(v1, v2, elementIndex);
    }

    template<class Derived>
    void LValueMatrix<Derived>::majorReduce(size_t v1, size_t v2, const ScalarType& factor) {
        if constexpr (MatrixOption::isColMatrix<Derived>()) {
            auto col1 = col(v1);
            col1 -= col(v2) * factor;
        }
        else {
            auto row1 = row(v1);
            row1 -= row(v2) * factor;
        }
    }

    template<class Derived>
    void LValueMatrix<Derived>::majorMulScalar(size_t v, const ScalarType& factor) {
        if constexpr (MatrixOption::isColMatrix<Derived>()) {
            auto c = col(v);
            c *= factor;
        }
        else {
            auto r = row(v);
            r *= factor;
        }
    }

    template<class Derived>
    void LValueMatrix<Derived>::majorSwap(size_t v1, size_t v2) {
        if constexpr (MatrixOption::isColMatrix<Derived>())
            Base::getDerived().swap_col(v1, v2);
        else
            Base::getDerived().swap_row(v1, v2);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::flatten(this auto&& self) {
        using Self = decltype(self);
        return FlattenL<Self>(std::forward<Self>(self));
    }

    template<class Derived>
    void LValueMatrix<Derived>::zeros() noexcept {
        const size_t maxMajor = Base::getMaxMajor();
        const size_t maxMinor = Base::getMaxMinor();
        for (size_t major = 0; major < maxMajor; ++major)
            for (size_t minor = 0; minor < maxMinor; ++minor)
                refFromMajorMinor(major, minor) = T(0);
    }

    template<class Derived>
    void LValueMatrix<Derived>::clamp_min(Tv minimum) {
        Base::getDerived().flatten().clamp_min(minimum);
    }

    template<class Derived>
    void LValueMatrix<Derived>::clamp_max(Tv maximum) {
        Base::getDerived().flatten().clamp_max(maximum);
    }

    template<class Derived>
    template<RNG R>
    void LValueMatrix<Derived>::random_uniform() {
        const size_t maxMajor = Base::getMaxMajor();
        const size_t maxMinor = Base::getMaxMinor();
        for (size_t major = 0; major < maxMajor; ++major)
            for (size_t minor = 0; minor < maxMinor; ++minor)
                refFromMajorMinor(major, minor) = ScalarType::template random_uniform<R>();
    }

    template<class Derived>
    template<RNG R>
    void LValueMatrix<Derived>::random_normal() {
        const size_t maxMajor = Base::getMaxMajor();
        const size_t maxMinor = Base::getMaxMinor();
        for (size_t major = 0; major < maxMajor; ++major)
            for (size_t minor = 0; minor < maxMinor; ++minor)
                refFromMajorMinor(major, minor) = ScalarType::template random_normal<R>();
    }

    template<class Derived>
    template<RNG R>
    void LValueMatrix<Derived>::random_any(auto& distribution) {
        const size_t maxMajor = Base::getMaxMajor();
        const size_t maxMinor = Base::getMaxMinor();
        for (size_t major = 0; major < maxMajor; ++major)
            for (size_t minor = 0; minor < maxMinor; ++minor)
                refFromMajorMinor(major, minor) = ScalarType::template random_any<R>(distribution);
    }

    template<class Derived>
    void LValueMatrix<Derived>::toIdentity() {
        assert(Base::getRow() == Base::getCol());
        const size_t order = Base::getRow();
        for (size_t i = 0; i < order; ++i)
            for (size_t j = 0; j < order; ++j)
                refFromMajorMinor(i, j) = Trv(i == j ? 1 : 0);
    }
    /**
     * Transforms like D^{-1}AD, \returns diagonal part of matrix D
     * Using exact minimum instead of the one from [1] is faster and usually works good, though it introduces some rounding error.
     *
     * Reference:
     * [1] arxiv:1401.5766; https://doi.org/10.48550/arxiv.1401.5766
     */
    template<class Derived>
    auto LValueMatrix<Derived>::balance() -> VectorND<T> {
        assert_balance();

        const size_t order = Base::getRow();
        auto& m = Base::getDerived();
        VectorND<T> result(order, 1);
        size_t numIteration = 0;
        while (true) {
            bool converged = true;
            for (size_t i = 0; i < order; ++i) {
                Trv normR = m.row(i).squaredNorm();
                Trv normC = m.col(i).squaredNorm();
                Trv fac = sqrt(normC / normR);
                Trv rep = reciprocal(fac);
                if ((fac + rep) >= Trv(2 / 0.95)) {
                    m.col(i) *= rep;
                    m.row(i) *= fac;
                    result[i] *= fac;
                    converged = false;
                }
            }

            if (converged || numIteration >= 32) // 32 usually large enough
                break;
            numIteration += 1;
        }
        return result;
    }

    template<class Derived>
    template<int GradOrder>
    auto LValueMatrix<Derived>::grads() const noexcept {
        return Base::template grads_impl<GradOrder>();
    }

    template<class Derived>
    auto LValueMatrix<Derived>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        assert(row < self.getRow() && col < self.getCol());
        return self.getDerived().data_ptr(row, col);
    }

    template<class Derived>
    decltype(auto) LValueMatrix<Derived>::refFromMajorMinor(this auto&& self, size_t major, size_t minor) noexcept {
        const size_t r = MatrixOption::rowFromMajorMinor<Derived>(major, minor);
        const size_t c = MatrixOption::colFromMajorMinor<Derived>(major, minor);
        assert(r < self.getRow() && c < self.getCol());
        return self[r, c];
    }

    template<class Derived>
    void LValueMatrix<Derived>::assert_balance() const noexcept {
        static_assert(!MatrixOption::isSymmMatrix<Derived>(), "[Error]: Unnecesary balancing on symm matrix");
        assert(Base::isSquare() && "[Error]: balance() requires a square matrix");
    }
}
