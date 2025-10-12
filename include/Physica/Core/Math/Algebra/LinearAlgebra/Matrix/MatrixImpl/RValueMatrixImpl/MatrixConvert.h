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

#include "../RValueMatrix.h"

namespace Physica {
    template<class M>
    class RealMatrix : public RValueMatrix<RealMatrix<M>> {
        using This = RealMatrix<M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        const M& mat;
    public:
        RealMatrix(const M& mat_) : mat(mat_) {}
        RealMatrix(const This&) = default;
        RealMatrix(This&&) noexcept = default;
        ~RealMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const { return mat.calc(row, col).real(); }
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const { return calc(row, col).value(); }
        /* Getters */
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
    };

    template<class M>
    class ImagMatrix : public RValueMatrix<ImagMatrix<M>> {
        using This = ImagMatrix<M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        const M& mat;
    public:
        ImagMatrix(const M& mat_) : mat(mat_) {}
        ImagMatrix(const This&) = default;
        ImagMatrix(This&&) noexcept = default;
        ~ImagMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const { return mat.calc(row, col).imag(); }
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const { return calc(row, col).value(); }
        /* Getters */
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
    };

    template<class M>
    class SquaredNormMatrix : public RValueMatrix<SquaredNormMatrix<M>> {
        using This = SquaredNormMatrix<M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::Trv;
    private:
        const M& mat;
    public:
        SquaredNormMatrix(const M& mat_) : mat(mat_) {}
        SquaredNormMatrix(const This&) = default;
        SquaredNormMatrix(This&&) noexcept = default;
        ~SquaredNormMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const { return mat.calc(row, col).squaredNorm(); }
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const { return mat.calc(row, col).value().squaredNorm(); }

        void reverse(const auto& grad) const noexcept;
        using Base::reverse;
        /* Getters */
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
    };

    template<class M>
    void SquaredNormMatrix<M>::reverse(const auto& grad) const noexcept {
        using U = decltype(grad);
        if constexpr (Scalar<U>)
            mat.reverse((Trv(2) * grad) * mat.values());
        else {
            static_assert(Matrix<U>);
            mat.reverse(hadamard(Trv(2) * mat.values(), grad));
        }
    }

    template<class M>
    class NormMatrix : public RValueMatrix<NormMatrix<M>> {
        using This = NormMatrix<M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        const M& mat;
    public:
        NormMatrix(const M& mat_) : mat(mat_) {}
        NormMatrix(const This&) = default;
        NormMatrix(This&&) noexcept = default;
        ~NormMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const { return mat.calc(row, col).norm(); }
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const { return mat.calc(row, col).value().norm(); }
        /* Getters */
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
    };

    template<class M>
    class ValueMatrix : public RValueMatrix<ValueMatrix<M>> {
        using This = ValueMatrix<M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
    private:
        const M& mat;
    public:
        ValueMatrix(const M& mat_) : mat(mat_) {}
        ValueMatrix(const This&) = default;
        ValueMatrix(This&&) noexcept = default;
        ~ValueMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const { return mat.calc_value(row, col); }
        [[nodiscard]] T calc_value(size_t row, size_t col) const { return calc(row, col); }
        /* Getters */
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
    };

    template<class M, int GradOrder>
    class GradMatrix : public RValueMatrix<GradMatrix<M, GradOrder>> {
        using This = GradMatrix<M, GradOrder>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        const M& mat;
    public:
        GradMatrix(const M& mat_) : mat(mat_) {}
        GradMatrix(const This&) = default;
        GradMatrix(This&&) noexcept = default;
        ~GradMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const { return mat.calc(row, col).template grad<GradOrder>(); }
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const { return calc(row, col).value(); }
        /* Getters */
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
    };
}

namespace Physica {
    template<class M>
    class Traits<RealMatrix<M>> {
    public:
        using ScalarType = M::ScalarType::RealType;
        constexpr static int Option = M::Option;
        constexpr static size_t RowAtCompile = M::RowAtCompile;
        constexpr static size_t ColAtCompile = M::ColAtCompile;
        constexpr static size_t SizeAtCompile = M::SizeAtCompile;
    };

    template<class M>
    class Traits<ImagMatrix<M>> : public Traits<RealMatrix<M>> {};

    template<class M>
    class Traits<SquaredNormMatrix<M>> : public Traits<RealMatrix<M>> {};

    template<class M>
    class Traits<NormMatrix<M>> : public Traits<RealMatrix<M>> {};

    template<class M>
    class Traits<ValueMatrix<M>> {
    public:
        using ScalarType = M::ScalarType::ValueType;
        constexpr static int Option = M::Option;
        constexpr static size_t RowAtCompile = M::RowAtCompile;
        constexpr static size_t ColAtCompile = M::ColAtCompile;
        constexpr static size_t SizeAtCompile = M::SizeAtCompile;
    };

    template<class M, int GradOrder>
    class Traits<GradMatrix<M, GradOrder>> {
        static_assert(M::ScalarType::isDiffable, "[Error]: Unnecessary toValueVector() call or toGradVector() call");
    public:
        using ScalarType = Internal::GradTypeHelper<typename M::ScalarType, GradOrder>::Type;
        constexpr static int Option = M::Option;
        constexpr static size_t RowAtCompile = M::RowAtCompile;
        constexpr static size_t ColAtCompile = M::ColAtCompile;
        constexpr static size_t SizeAtCompile = M::SizeAtCompile;
    };
}
