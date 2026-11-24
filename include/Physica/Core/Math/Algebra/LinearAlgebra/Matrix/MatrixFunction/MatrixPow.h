/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RValueMatrix.h"

namespace Physica {
    template<Matrix M, Vector V> class MatPowVecProd;

    template<Matrix M>
    class MatrixPow : public RValueMatrix<MatrixPow<M>> {
        static_assert(std::is_reference_v<M>);

        using This = MatrixPow<M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<M> m;
        int power;
    public:
        MatrixPow(M m_, int power_);
        MatrixPow(const This&) = default;
        MatrixPow(This&&) noexcept = default;
        ~MatrixPow() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;

        template<Vector V>
        [[nodiscard]] auto operator*(this auto&&, V&& v) noexcept;
        /* Operations */
        [[nodiscard]] T calc(size_t, size_t) const { noImpl("MatrixPow::calc() is low performance and should be avoided"); }

        [[nodiscard]] auto transpose() const noexcept;
        [[nodiscard]] auto hermite() const noexcept;
        /* Getters */
        [[nodiscard]] const auto& getMatrix() const noexcept { return m; }
        [[nodiscard]] constexpr size_t getOrder() const noexcept { return m.getOrder(); }
        [[nodiscard]] constexpr int getPower() const noexcept { return power; }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isFastAssign() noexcept { return true; }
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept { return true; }
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return MatrixMajor::BothMajor; }
    };

    template<Matrix M>
    MatrixPow<M>::MatrixPow(M m_, int power_) : m(std::forward<M>(m_)), power(power_) {
        assert(m.isSquare());
    }

    template<Matrix M>
    auto MatrixPow<M>::transpose() const noexcept {
        return pow(m.transpose(), power);
    }

    template<Matrix M>
    auto MatrixPow<M>::hermite() const noexcept {
        return pow(m.hermite(), power);
    }

    template<Matrix M>
    template<Vector V>
    auto MatrixPow<M>::operator*(this auto&& self, V&& v) noexcept {
        using Self = decltype(self);
        return MatPowVecProd<Self, V&&>(std::forward<Self>(self), std::forward<V>(v));
    }

    template<Matrix M>
    [[nodiscard]] auto pow(M&& m, int power) noexcept {
        return MatrixPow<M&&>(std::forward<M>(m), power);
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<MatrixPow<M>> : public Traits<std::remove_cvref_t<M>> {};
}

#include "MatPowVecProd.h"
