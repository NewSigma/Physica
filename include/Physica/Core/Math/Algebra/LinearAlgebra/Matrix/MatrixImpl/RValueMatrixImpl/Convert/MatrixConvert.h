/*
 * Copyright 2022-2026 Weibo He.
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

#include "../../RValueMatrix.h"

namespace Physica {
    template<class M>
    class RealMatrix : public RValueMatrix<RealMatrix<M>> {
        using This = RealMatrix<M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        decay_rvalue_t<M> mat;
    public:
        explicit RealMatrix(M&& mat_) : mat(std::forward<M>(mat_)) {}
        RealMatrix(const This&) = default;
        RealMatrix(This&&) noexcept = default;
        ~RealMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const { return mat.calc(row, col).real(); }

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
        [[nodiscard]] size_t getOrder() const { return mat.getOrder(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept { return std::remove_cvref_t<M>::isStaticSquare(); }
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept { return std::remove_cvref_t<M>::getRowAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept { return std::remove_cvref_t<M>::getColAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return std::remove_cvref_t<M>::getMajor(); }
    };

    template<class M>
    auto RealMatrix<M>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat).values().reals();
    }

    template<class M>
    class ImagMatrix : public RValueMatrix<ImagMatrix<M>> {
        using This = ImagMatrix<M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        decay_rvalue_t<M> mat;
    public:
        explicit ImagMatrix(M&& mat_) : mat(std::forward<M>(mat_)) {}
        ImagMatrix(const This&) = default;
        ImagMatrix(This&&) noexcept = default;
        ~ImagMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const { return mat.calc(row, col).imag(); }

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
        [[nodiscard]] size_t getOrder() const { return mat.getOrder(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept { return std::remove_cvref_t<M>::isStaticSquare(); }
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept { return std::remove_cvref_t<M>::getRowAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept { return std::remove_cvref_t<M>::getColAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return std::remove_cvref_t<M>::getMajor(); }
    };

    template<class M>
    auto ImagMatrix<M>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat).values().imags();
    }

    template<class M>
    class SquaredNormMatrix : public RValueMatrix<SquaredNormMatrix<M>> {
        using This = SquaredNormMatrix<M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::Trv;
    private:
        decay_rvalue_t<M> mat;
    public:
        explicit SquaredNormMatrix(M&& mat_) : mat(std::forward<M>(mat_)) {}
        SquaredNormMatrix(const This&) = default;
        SquaredNormMatrix(This&&) noexcept = default;
        ~SquaredNormMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const { return mat.calc(row, col).squaredNorm(); }

        void reverse(const auto& grad) const noexcept;
        using Base::reverse;

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
        [[nodiscard]] size_t getOrder() const { return mat.getOrder(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept { return std::remove_cvref_t<M>::isStaticSquare(); }
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept { return std::remove_cvref_t<M>::getRowAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept { return std::remove_cvref_t<M>::getColAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return std::remove_cvref_t<M>::getMajor(); }
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
    auto SquaredNormMatrix<M>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat).values().squaredNorms();
    }

    template<class M>
    class NormMatrix : public RValueMatrix<NormMatrix<M>> {
        using This = NormMatrix<M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        decay_rvalue_t<M> mat;
    public:
        explicit NormMatrix(M&& mat_) : mat(std::forward<M>(mat_)) {}
        NormMatrix(const This&) = default;
        NormMatrix(This&&) noexcept = default;
        ~NormMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] CoDiff<T> calc(size_t row, size_t col) const { return mat.calc(row, col).norm(); }

        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
        [[nodiscard]] size_t getOrder() const { return mat.getOrder(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept { return std::remove_cvref_t<M>::isStaticSquare(); }
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept { return std::remove_cvref_t<M>::getRowAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept { return std::remove_cvref_t<M>::getColAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return std::remove_cvref_t<M>::getMajor(); }
    };

    template<class M>
    auto NormMatrix<M>::values(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat).values().norms();
    }

    template<class M>
    class ValueMatrix : public RValueMatrix<ValueMatrix<M>> {
        using This = ValueMatrix<M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<M> mat;
    public:
        explicit ValueMatrix(M&& mat_) : mat(std::forward<M>(mat_)) {}
        ValueMatrix(const This&) = default;
        ValueMatrix(This&&) noexcept = default;
        ~ValueMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const { return mat.calc(row, col).value(); }
        [[nodiscard]] T calc_value(size_t row, size_t col) const { return calc(row, col); }
        /* Getters */
        [[nodiscard]] size_t getRow() const { return mat.getRow(); }
        [[nodiscard]] size_t getCol() const { return mat.getCol(); }
        [[nodiscard]] size_t getOrder() const { return mat.getOrder(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept { return std::remove_cvref_t<M>::isStaticSquare(); }
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept { return std::remove_cvref_t<M>::getRowAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept { return std::remove_cvref_t<M>::getColAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return std::remove_cvref_t<M>::getMajor(); }
    };

    template<class M, int GradOrder>
    class GradMatrix : public RValueMatrix<GradMatrix<M, GradOrder>> {
        using This = GradMatrix<M, GradOrder>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        decay_rvalue_t<M> mat;
    public:
        explicit GradMatrix(M&& mat_) : mat(std::forward<M>(mat_)) {}
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
        [[nodiscard]] size_t getOrder() const { return mat.getOrder(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept { return std::remove_cvref_t<M>::isStaticSquare(); }
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept { return std::remove_cvref_t<M>::getRowAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept { return std::remove_cvref_t<M>::getColAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return std::remove_cvref_t<M>::getMajor(); }
    };
}

namespace Physica {
    template<class M>
    class Traits<RealMatrix<M>> {
        using M1 = std::remove_cvref_t<M>;
    public:
        using ScalarType = typename M1::ScalarType::RealType;
    };

    template<class M>
    class Traits<ImagMatrix<M>> : public Traits<RealMatrix<M>> {};

    template<class M>
    class Traits<SquaredNormMatrix<M>> : public Traits<RealMatrix<M>> {};

    template<class M>
    class Traits<NormMatrix<M>> : public Traits<RealMatrix<M>> {};

    template<class M>
    class Traits<ValueMatrix<M>> {
        using M1 = std::remove_cvref_t<M>;
    public:
        using ScalarType = typename M1::ScalarType::ValueType;
    };

    template<class M, int GradOrder>
    class Traits<GradMatrix<M, GradOrder>> {
        using M1 = std::remove_cvref_t<M>;
        static_assert(M1::isDiffable(), "[Error]: Redundant GradMatrix");
    public:
        using ScalarType = M1::ScalarType::template GradWithOrder<GradOrder>;
    };
}
