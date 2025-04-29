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

#include "../LValueMatrix.h"

namespace Physica {
    template<Vector T, int MatrixMajor, size_t Row, size_t Col>
    class RValueReshapedVector;

    template<Vector T, int MatrixMajor, size_t Row, size_t Col>
    class LValueReshapedVector : public LValueMatrix<LValueReshapedVector<T, MatrixMajor, Row, Col>> {
        using This = LValueReshapedVector<T, MatrixMajor, Row, Col>;
        using Base = LValueMatrix<This>;
    public:
        using typename Base::ScalarType;
    protected:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    private:
        T& v;
        size_t r;
        size_t c;
    public:
        LValueReshapedVector(T& v_, size_t r_, size_t c_);
        LValueReshapedVector(const This&) = default;
        LValueReshapedVector(This&&) noexcept = default;
        ~LValueReshapedVector() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        using Base::operator=;
        /* Operations */
        void resize(size_t row, size_t col);
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept;
        [[nodiscard]] size_t getCol() const noexcept;
        [[nodiscard]] inline PtrTy data_ptr(size_t row, size_t col) noexcept;
        [[nodiscard]] inline ConstPtrTy data_ptr(size_t row, size_t col) const noexcept;
        [[nodiscard]] ScalarType sum() const { return v.sum(); }
    };

    template<Vector T, int MatrixMajor, size_t Row, size_t Col>
    LValueReshapedVector<T, MatrixMajor, Row, Col>::LValueReshapedVector(T& v_, size_t r_, size_t c_)
            : v(v_), r(r_), c(c_) {
        assert(r == Row || Row == Dynamic);
        assert(c == Col || Col == Dynamic);
        assert(r * c == v.getLength());
    }

    template<Vector T, int MatrixMajor, size_t Row, size_t Col>
    void LValueReshapedVector<T, MatrixMajor, Row, Col>::resize(
            [[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
        assert(row == getRow());
        assert(col == getCol());
    }

    template<Vector T, int MatrixMajor, size_t Row, size_t Col>
    size_t LValueReshapedVector<T, MatrixMajor, Row, Col>::getRow() const noexcept {
        if constexpr (Row != Dynamic)
            return Row;
        return r;
    }

    template<Vector T, int MatrixMajor, size_t Row, size_t Col>
    size_t LValueReshapedVector<T, MatrixMajor, Row, Col>::getCol() const noexcept {
        if constexpr (Col != Dynamic)
            return Col;
        return c;
    }

    template<Vector T, int MatrixMajor, size_t Row, size_t Col>
    inline auto LValueReshapedVector<T, MatrixMajor, Row, Col>::data_ptr(size_t row, size_t col) noexcept -> PtrTy {
        assert(row < getRow() && col < getCol());
        if constexpr (MatrixOption::isColMatrix<This>())
            return v.data_ptr(col * getRow() + row);
        else
            return v.data_ptr(row * getCol() + col);
    }

    template<Vector T, int MatrixMajor, size_t Row, size_t Col>
    inline auto LValueReshapedVector<T, MatrixMajor, Row, Col>::data_ptr(size_t row, size_t col) const noexcept -> ConstPtrTy {
        return const_cast<This&>(*this).data_ptr(row, col);
    }

    template<class Derived>
    template<Matrix M>
    auto LValueVector<Derived>::reshape(const M& mat) noexcept {
        using ResultType = LValueReshapedVector<Derived, MatrixOption::getMajor<M>(), M::RowAtCompile, M::ColAtCompile>;
        return ResultType(Base::getDerived(), mat.getRow(), mat.getCol());
    }

    template<class Derived>
    template<Matrix M>
    const auto LValueVector<Derived>::reshape(const M& mat) const noexcept {
        using ResultType = LValueReshapedVector<Derived, MatrixOption::getMajor<M>(), M::RowAtCompile, M::ColAtCompile>;
        return ResultType(Base::getConstCastDerived(), mat.getRow(), mat.getCol());
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto LValueVector<Derived>::reshape_col(size_t row, size_t col) noexcept {
        return LValueReshapedVector<Derived, MatrixOption::Col, Row, Col>(Base::getDerived(), row, col);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    const auto LValueVector<Derived>::reshape_col(size_t row, size_t col) const noexcept {
        return LValueReshapedVector<Derived, MatrixOption::Col, Row, Col>(Base::getConstCastDerived(), row, col);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto LValueVector<Derived>::reshape_row(size_t row, size_t col) noexcept {
        return LValueReshapedVector<Derived, MatrixOption::Row, Row, Col>(Base::getDerived(), row, col);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    const auto LValueVector<Derived>::reshape_row(size_t row, size_t col) const noexcept {
        return LValueReshapedVector<Derived, MatrixOption::Row, Row, Col>(Base::getConstCastDerived(), row, col);
    }
}

namespace Physica {
    template<Vector T, int MatrixMajor, size_t Row, size_t Col>
    class Traits<LValueReshapedVector<T, MatrixMajor, Row, Col>>
            : public Traits<RValueReshapedVector<T, MatrixMajor, Row, Col>> {};
}
