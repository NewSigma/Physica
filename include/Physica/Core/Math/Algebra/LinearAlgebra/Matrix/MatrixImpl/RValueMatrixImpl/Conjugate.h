/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "../RValueMatrix.h"

namespace Physica {
    template<Matrix M>
    class Conjugate<M> : public RValueMatrix<Conjugate<M>> {
        using This = Conjugate<M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<M> m;
    public:
        Conjugate(M&& m);
        Conjugate(const This&) = default;
        Conjugate(This&&) noexcept = default;
        ~Conjugate() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const { return m.calc(row, col).conjugate(); }

        [[nodiscard]] decltype(auto) transpose(this auto&&) noexcept;
        [[nodiscard]] auto&& conjugate(this auto&&) noexcept;
        [[nodiscard]] decltype(auto) hermite(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] constexpr size_t getRow() const noexcept { return m.getRow(); }
        [[nodiscard]] constexpr size_t getCol() const noexcept { return m.getCol(); }
        [[nodiscard]] constexpr size_t getOrder() const noexcept { return m.getOrder(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSymm() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isStaticHermite() noexcept;
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept;
    };

    template<Matrix M>
    Conjugate<M>::Conjugate(M&& m) : m(std::forward<M>(m)) {
        static_assert(m.isComplex(), "[Error]: Unnecessary conjugate call on real matrix");
    }

    template<Matrix M>
    decltype(auto) Conjugate<M>::transpose(this auto&& self) noexcept {
        using Self = decltype(self);
        return propagate_rvalue_reference<Self, M>(self.m).hermite();
    }

    template<Matrix M>
    auto&& Conjugate<M>::conjugate(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.m);
    }

    template<Matrix M>
    decltype(auto) Conjugate<M>::hermite(this auto&& self) noexcept {
        using Self = decltype(self);
        return propagate_rvalue_reference<Self, M>(self.m).transpose();
    }

    template<Matrix M>
    __host__ __device__ consteval bool Conjugate<M>::isStaticSymm() noexcept {
        return std::remove_cvref_t<M>::isStaticSymm();
    }

    template<Matrix M>
    __host__ __device__ consteval bool Conjugate<M>::isStaticHermite() noexcept {
        return std::remove_cvref_t<M>::isStaticHermite();
    }

    template<Matrix M>
    __host__ __device__ consteval bool Conjugate<M>::isStaticSquare() noexcept {
        return std::remove_cvref_t<M>::isStaticSquare();
    }

    template<Matrix M>
    __host__ __device__ consteval size_t Conjugate<M>::getRowAtCompile() noexcept {
        return std::remove_cvref_t<M>::getRowAtCompile();
    }

    template<Matrix M>
    __host__ __device__ consteval size_t Conjugate<M>::getColAtCompile() noexcept {
        return std::remove_cvref_t<M>::getColAtCompile();
    }

    template<Matrix M>
    __host__ __device__ consteval int Conjugate<M>::getMajor() noexcept {
        return std::remove_cvref_t<M>::getMajor();
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<Conjugate<M>> : public Traits<M> {};
}
