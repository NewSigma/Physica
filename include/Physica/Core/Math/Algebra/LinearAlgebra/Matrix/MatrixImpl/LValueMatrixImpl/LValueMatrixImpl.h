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
    template<Scalar T>
    Derived& LValueMatrix<Derived>::operator=(const T& x) requires(!isReverseDiff || !ReverseDiff<T>) {
        ScalarType::template assert_assign<T>();
        const size_t maxMajor = Base::getMaxMajor();
        const size_t maxMinor = Base::getMaxMinor();
        for (size_t i = 0; i < maxMajor; ++i) {
            for (size_t j = 0; j < maxMinor; ++j) {
                if constexpr (ReverseDiff<T>)
                    refFromMajorMinor(i, j) = x.value();
                else
                    refFromMajorMinor(i, j) = x;
            }
        }
        return Base::getDerived();
    }

    template<class Derived>
    void LValueMatrix<Derived>::operator+=(const Scalar auto& x) {
        auto& m = Base::getDerived();
        (m + x).assign(m);
    }

    template<class Derived>
    void LValueMatrix<Derived>::operator-=(const Scalar auto& x) {
        auto& m = Base::getDerived();
        (m - x).assign(m);
    }

    template<class Derived>
    void LValueMatrix<Derived>::operator*=(const Scalar auto& x) {
        auto& m = Base::getDerived();
        (m * x).assign(m);
    }

    template<class Derived>
    void LValueMatrix<Derived>::operator/=(const Scalar auto& x) {
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
    auto LValueMatrix<Derived>::row(size_t r) noexcept {
        return RowVector(Base::getDerived(), r, 0, Base::getCol());
    }

    template<class Derived>
    const auto LValueMatrix<Derived>::row(size_t r) const noexcept {
        return RowVector(Base::getConstCastDerived(), r, 0, Base::getCol());
    }

    template<class Derived>
    auto LValueMatrix<Derived>::col(size_t c) noexcept {
        return ColVector(Base::getDerived(), 0, Base::getRow(), c);
    }

    template<class Derived>
    const auto LValueMatrix<Derived>::col(size_t c) const noexcept {
        return ColVector(Base::getConstCastDerived(), 0, Base::getRow(), c);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::rows(size_t fromRow, size_t rowCount) noexcept {
        return BlockType(Base::getDerived(), fromRow, rowCount, 0, Base::getCol());
    }

    template<class Derived>
    const auto LValueMatrix<Derived>::rows(size_t fromRow, size_t rowCount) const noexcept {
        return BlockType(Base::getConstCastDerived(), fromRow, rowCount, 0, Base::getCol());
    }

    template<class Derived>
    auto LValueMatrix<Derived>::topRows(size_t to) noexcept {
        return BlockType(Base::getDerived(), 0, to, 0, Base::getCol());
    }

    template<class Derived>
    const auto LValueMatrix<Derived>::topRows(size_t to) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, to, 0, Base::getCol());
    }

    template<class Derived>
    auto LValueMatrix<Derived>::bottomRows(size_t from) noexcept {
        return BlockType(Base::getDerived(), from, Base::getRow() - from, 0, Base::getCol());
    }

    template<class Derived>
    const auto LValueMatrix<Derived>::bottomRows(size_t from) const noexcept {
        return BlockType(Base::getConstCastDerived(), from, Base::getRow() - from, 0, Base::getCol());
    }

    template<class Derived>
    auto LValueMatrix<Derived>::cols(size_t fromCol, size_t colCount) noexcept {
        return BlockType(Base::getDerived(), 0, Base::getRow(), fromCol, colCount);
    }

    template<class Derived>
    const auto LValueMatrix<Derived>::cols(size_t fromCol, size_t colCount) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, Base::getRow(), fromCol, colCount);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::leftCols(size_t to) noexcept {
        return BlockType(Base::getDerived(), 0, Base::getRow(), 0, to);
    }

    template<class Derived>
    const auto LValueMatrix<Derived>::leftCols(size_t to) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, Base::getRow(), 0, to);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::rightCols(size_t from) noexcept {
        return BlockType(Base::getDerived(), 0, Base::getRow(), from, Base::getCol() - from);
    }

    template<class Derived>
    const auto LValueMatrix<Derived>::rightCols(size_t from) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, Base::getRow(), from, Base::getCol() - from);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::topLeftCorner(size_t toRow, size_t toCol) noexcept {
        return BlockType(Base::getDerived(), 0, toRow, 0, toCol);
    }

    template<class Derived>
    const auto LValueMatrix<Derived>::topLeftCorner(size_t toRow, size_t toCol) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, toRow, 0, toCol);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::topLeftCorner(size_t to) noexcept {
        return BlockType(Base::getDerived(), 0, to, 0, to);
    }

    template<class Derived>
    const auto LValueMatrix<Derived>::topLeftCorner(size_t to) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, to, 0, to);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::topRightCorner(size_t toRow, size_t fromCol) noexcept {
        return BlockType(Base::getDerived(), 0, toRow, fromCol, Base::getRow() - fromCol);
    }

    template<class Derived>
    const auto LValueMatrix<Derived>::topRightCorner(size_t toRow, size_t fromCol) const noexcept {
        return BlockType(Base::getConstCastDerived(), 0, toRow, fromCol, Base::getRow() - fromCol);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::bottomLeftCorner(size_t fromRow, size_t toCol) noexcept {
        return BlockType(Base::getDerived(), fromRow, Base::getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    const auto LValueMatrix<Derived>::bottomLeftCorner(size_t fromRow, size_t toCol) const noexcept {
        return BlockType(Base::getConstCastDerived(), fromRow, Base::getRow() - fromRow, 0, toCol);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::bottomRightCorner(size_t fromRow, size_t fromCol) noexcept {
        return BlockType(Base::getDerived(), fromRow, Base::getRow() - fromRow, fromCol, Base::getCol() - fromCol);
    }

    template<class Derived>
    const auto LValueMatrix<Derived>::bottomRightCorner(size_t fromRow, size_t fromCol) const noexcept {
        return BlockType(Base::getConstCastDerived(), fromRow, Base::getRow() - fromRow, fromCol, Base::getCol() - fromCol);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::bottomRightCorner(size_t from) noexcept {
        return BlockType(Base::getDerived(), from, Base::getRow() - from, from, Base::getCol() - from);
    }

    template<class Derived>
    const auto LValueMatrix<Derived>::bottomRightCorner(size_t from) const noexcept {
        return BlockType(Base::getConstCastDerived(), from, Base::getRow() - from, from, Base::getCol() - from);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        return BlockType(Base::getDerived(), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    const auto LValueMatrix<Derived>::block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const noexcept {
        return BlockType(Base::getConstCastDerived(), fromRow, rowCount, fromCol, colCount);
    }

    template<class Derived>
    auto LValueMatrix<Derived>::diag() noexcept {
        assert(Base::isSquare());
        return DiagVector<Derived, true>(Base::getDerived());
    }

    template<class Derived>
    const auto LValueMatrix<Derived>::diag() const noexcept {
        assert(Base::isSquare());
        return DiagVector<Derived, true>(Base::getConstCastDerived());
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
        auto row2 = row(r2);
        row2 -= row(r1) * factor;
        matrix[r2, elementIndex] = 0;
    }

    template<class Derived>
    void LValueMatrix<Derived>::colReduce(size_t c1, size_t c2, size_t elementIndex) {
        Derived& matrix = Base::getDerived();
        assert(abs(matrix(elementIndex, c1)) > ScalarType(std::numeric_limits<ScalarType>::min()));
        const ScalarType factor = matrix(c2, elementIndex) / matrix(c1, elementIndex);
        auto col2 = col(c2);
        col2 -= col(c1) * factor;
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
    auto LValueMatrix<Derived>::flatten() {
        return FlattenL<Derived>(Base::getDerived());
    }

    template<class Derived>
    const auto LValueMatrix<Derived>::flatten() const {
        return FlattenL<Derived>(Base::getConstCastDerived());
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
