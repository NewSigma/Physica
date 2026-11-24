/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "../RValueMatrix.h"

namespace Physica {
    template<class T>
    struct remove_transpose {
        using Type = T;
    };

    template<class T>
    struct remove_transpose<Transpose<T>> {
        using Type = T;
    };

    template<Matrix M>
    class Transpose<M> : public RValueMatrix<Transpose<M>> {
        using This = Transpose<M>;
        using Base = RValueMatrix<This>;
    public:
        using typename Base::T;
        using typename Base::Tv;
    private:
        decay_rvalue_t<M> mat;
    public:
        constexpr Transpose(M&& mat_) : mat(std::forward<M>(mat_)) {}
        Transpose(const This&) = default;
        Transpose(This&&) noexcept = default;
        ~Transpose() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] constexpr T calc(size_t row, size_t col) const;

        void reverse(const Matrix auto& grad) const noexcept;
        using Base::reverse;

        [[nodiscard]] constexpr auto&& transpose(this auto&&) noexcept;
        [[nodiscard]] decltype(auto) conjugate(this auto&&) noexcept;
        [[nodiscard]] decltype(auto) hermite(this auto&&) noexcept;
        [[nodiscard]] decltype(auto) flatten(this auto&&) noexcept;
        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] constexpr size_t getRow() const noexcept { return mat.getCol(); }
        [[nodiscard]] constexpr size_t getCol() const noexcept { return mat.getRow(); }
        [[nodiscard]] constexpr size_t getOrder() const noexcept { return mat.getOrder(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSymm() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isStaticHermite() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getMinorStrideAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept;
    };

    template<Matrix M>
    constexpr auto Transpose<M>::calc(size_t row, size_t col) const -> T {
        return mat.calc(col, row);
    }

    template<Matrix M>
    void Transpose<M>::reverse(const Matrix auto& grad) const noexcept {
        static_assert(Base::isReverseDiff());
        mat.reverse(grad.transpose());
    }

    template<Matrix M>
    constexpr auto&& Transpose<M>::transpose(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }

    template<Matrix M>
    decltype(auto) Transpose<M>::conjugate(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (T::isComplex())
            return propagate_rvalue_reference<Self, M>(self.mat).hermite();
        else
            return std::forward<Self>(self);
    }

    template<Matrix M>
    decltype(auto) Transpose<M>::hermite(this auto&& self) noexcept {
        using Self = decltype(self);
        return propagate_rvalue_reference<Self, M>(self.mat).conjugate();
    }

    template<Matrix M>
    decltype(auto) Transpose<M>::flatten(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).transpose().flatten();
    }

    template<Matrix M>
    auto Transpose<M>::values(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).transpose().values().transpose();
    }

    template<Matrix M>
    __host__ __device__ consteval bool Transpose<M>::isStaticSymm() noexcept {
        static_assert(!std::remove_cvref_t<M>::isStaticSymm(), "[Error]: Do not construct Transpose<> for a statically symmetric matrix");
        return false;
    }

    template<Matrix M>
    __host__ __device__ consteval bool Transpose<M>::isStaticHermite() noexcept {
        return std::remove_cvref_t<M>::isStaticHermite();
    }

    template<Matrix M>
    __host__ __device__ consteval bool Transpose<M>::isStaticSquare() noexcept {
        return std::remove_cvref_t<M>::isStaticSquare();
    }

    template<Matrix M>
    __host__ __device__ consteval size_t Transpose<M>::getRowAtCompile() noexcept {
        return std::remove_cvref_t<M>::getColAtCompile();
    }

    template<Matrix M>
    __host__ __device__ consteval size_t Transpose<M>::getColAtCompile() noexcept {
        return std::remove_cvref_t<M>::getRowAtCompile();
    }

    template<Matrix M>
    __host__ __device__ consteval size_t Transpose<M>::getMinorStrideAtCompile() noexcept {
        return std::remove_cvref_t<M>::getMinorStrideAtCompile();
    }

    template<Matrix M>
    __host__ __device__ consteval int Transpose<M>::getMajor() noexcept {
        constexpr static int OtherMajor = MatrixMajor::isColMatrix<M>() ? MatrixMajor::Row : MatrixMajor::Col;
        return MatrixMajor::isBothMajor<M>() ? MatrixMajor::BothMajor : OtherMajor;
    }

    template<Vector V>
    class Transpose<V> : public RValueMatrix<Transpose<V>> {
        using This = Transpose<V>;
        using Base = RValueMatrix<This>;
    public:
        using typename Base::T;
        using typename Base::Tv;
    private:
        decay_rvalue_t<V> vec;
    public:
        constexpr explicit Transpose(V&& vec_) : vec(std::forward<V>(vec_)) {}
        Transpose(const This&) = default;
        Transpose(This&&) noexcept = default;
        ~Transpose() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] constexpr T calc(size_t row, size_t col) const;

        [[nodiscard]] constexpr auto&& transpose(this auto&&) noexcept;
        [[nodiscard]] decltype(auto) conjugate(this auto&&) noexcept;
        [[nodiscard]] decltype(auto) hermite(this auto&&) noexcept;
        [[nodiscard]] decltype(auto) flatten(this auto&&) noexcept;
        [[nodiscard]] auto values(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] constexpr size_t getCol() const noexcept { return vec.getLength(); }
        [[nodiscard]] constexpr size_t getOrder() const noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return MatrixMajor::Row; }
    };

    template<Vector V>
    constexpr auto Transpose<V>::calc([[maybe_unused]] size_t row, size_t col) const -> T {
        assert(row == 0);
        return vec.calc(col);
    }

    template<Vector V>
    constexpr auto&& Transpose<V>::transpose(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.vec);
    }

    template<Vector V>
    decltype(auto) Transpose<V>::conjugate(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (T::isComplex())
            return propagate_rvalue_reference<Self, V>(self.vec).hermite();
        else
            return std::forward<Self>(self);
    }

    template<Vector V>
    decltype(auto) Transpose<V>::hermite(this auto&& self) noexcept {
        using Self = decltype(self);
        return propagate_rvalue_reference<Self, V>(self.vec).conjugate();
    }

    template<Vector V>
    decltype(auto) Transpose<V>::flatten(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).transpose();
    }

    template<Vector V>
    auto Transpose<V>::values(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).transpose().values().transpose();
    }

    template<Vector V>
    constexpr size_t Transpose<V>::getOrder() const noexcept {
        assert(Base::isSquare() && "[Error]: getOrder() assumes square matrix");
        return 1;
    }

    template<Vector V>
    __host__ __device__ consteval size_t Transpose<V>::getRowAtCompile() noexcept {
        return 1;
    }

    template<Vector V>
    __host__ __device__ consteval size_t Transpose<V>::getColAtCompile() noexcept {
        return std::remove_cvref_t<V>::getSizeAtCompile();
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<Transpose<M>> {
    public:
        using ScalarType = std::remove_cvref_t<M>::ScalarType;
    };

    template<Vector V>
    class Traits<Transpose<V>> {
    public:
        using ScalarType = std::remove_cvref_t<V>::ScalarType;
    };
}
