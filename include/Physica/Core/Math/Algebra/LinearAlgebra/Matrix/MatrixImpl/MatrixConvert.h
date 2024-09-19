/*
 * Copyright 2022-2024 Weibo He.
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

#include "RValueMatrix.h"

namespace Physica::Core {
    template<class MatrixType>
    class RealMatrix : public RValueMatrix<RealMatrix<MatrixType>> {
        using Base = RValueMatrix<RealMatrix<MatrixType>>;
        using This = RealMatrix<MatrixType>;
    public:
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
    public:
        RealMatrix(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}
        RealMatrix(const This&) = delete;
        RealMatrix(This&&) noexcept = delete;
        ~RealMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return mat.calc(row, col).getReal(); }
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    class ValueMatrix : public RValueMatrix<ValueMatrix<MatrixType>> {
        using Base = RValueMatrix<ValueMatrix<MatrixType>>;
        using This = ValueMatrix<MatrixType>;
    public:
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
    public:
        ValueMatrix(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}
        ValueMatrix(const This&) = delete;
        ValueMatrix(This&&) noexcept = delete;
        ~ValueMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return mat.calc(row, col).getValue(); }
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    class GradMatrix : public RValueMatrix<GradMatrix<MatrixType>> {
        using Base = RValueMatrix<GradMatrix<MatrixType>>;
        using This = GradMatrix<MatrixType>;
    public:
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
    public:
        GradMatrix(const RValueMatrix<MatrixType>& mat_) : mat(mat_.getDerived()) {}
        GradMatrix(const This&) = delete;
        GradMatrix(This&&) noexcept = delete;
        ~GradMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const { return mat.calc(row, col).getGrad(); }
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getColumn() const { return mat.getColumn(); }
    };

    template<class MatrixType>
    [[nodiscard]] auto toRealMatrix(const RValueMatrix<MatrixType>& m) {
        return RealMatrix(m);
    }

    template<class MatrixType>
    [[nodiscard]] auto toValueMatrix(const RValueMatrix<MatrixType>& m) {
        return ValueMatrix(m);
    }

    template<class MatrixType>
    [[nodiscard]] auto toGradMatrix(const RValueMatrix<MatrixType>& m) {
        return GradMatrix(m);
    }
}

namespace Physica {
    template <class MatrixType>
    class Traits<Core::RealMatrix<MatrixType>> {
    public:
        using ScalarType = typename MatrixType::ScalarType::RealType;
        constexpr static int Option = MatrixType::Option;
        constexpr static size_t RowAtCompile = MatrixType::RowAtCompile;
        constexpr static size_t ColumnAtCompile = MatrixType::ColumnAtCompile;
        constexpr static size_t SizeAtCompile = MatrixType::SizeAtCompile;
    };

    template <class MatrixType>
    class Traits<Core::ValueMatrix<MatrixType>> {
        using T = typename MatrixType::ScalarType;
        static_assert(T::isDifferentiable, "[Error]: Unnecessary toValueVector() call or toGradVector() call");
    public:
        using ScalarType = typename T::PlainScalar;
        constexpr static int Option = MatrixType::Option;
        constexpr static size_t RowAtCompile = MatrixType::RowAtCompile;
        constexpr static size_t ColumnAtCompile = MatrixType::ColumnAtCompile;
        constexpr static size_t SizeAtCompile = MatrixType::SizeAtCompile;
    };

    template <class MatrixType>
    class Traits<Core::GradMatrix<MatrixType>> : public Traits<Core::ValueMatrix<MatrixType>> {};
}
