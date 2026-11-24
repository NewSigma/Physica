/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "MatrixTrig.h"

namespace Physica {
    template<Matrix M> requires(instanceof_tx<M, MatrixTrig>)
    class Inverse<M> : public RValueMatrix<Inverse<M>> {
        using This = Inverse<M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
        using typename Base::Trv;
        using typename Base::Tc;
    private:
        decay_rvalue_t<M> trig;
    public:
        explicit Inverse(M&& trig);
        Inverse(const This&) = default;
        Inverse(This&&) noexcept = default;
        ~Inverse() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto& target) const;
        void assign_base(Matrix auto& target) const noexcept;
        void assign_mkl(Matrix auto& target) const;
        /* Getters */
        [[nodiscard]] auto&& getExpr(this auto&&) noexcept;
        [[nodiscard]] constexpr size_t getOrder() const noexcept { return trig.getOrder(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept { return true; }
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept;
    };

    template<Matrix M> requires(instanceof_tx<M, MatrixTrig>)
    Inverse<M>::Inverse(M&& trig) : trig(std::forward<M>(trig)) {
        assert(trig.isSquare() && "[Error]: inv() requires square matrix");
    }

    template<Matrix M> requires(instanceof_tx<M, MatrixTrig>)
    void Inverse<M>::assign(Matrix auto& target) const {
        using Expr = Traits<M>::ExprType;
        constexpr bool UseMKL = HasMKL() && Internal::EnableLAPACK<Expr, decltype(target)>::value && target.getMinorStrideAtCompile() == 1;
        if constexpr (UseMKL)
            assign_mkl(target);
        else
            assign_base(target);
    }

    template<Matrix M> requires(instanceof_tx<M, MatrixTrig>)
    void Inverse<M>::assign_base(Matrix auto& target) const noexcept {
        constexpr bool Unit = Traits<M>::Unit;
        target.assert_assign(*this);
        target.zeros();
        if constexpr (Unit)
            target.diag() = Trv(1);
        else
            target.diag() = reciprocal(trig.diag());

        if constexpr (Traits<M>::Upper) {
            for (size_t i = getOrder(); i > 1; --i) {
                const size_t k = i - 1;
                auto corner = target.topRightCorner(k, k);
                auto head = target.row(k).tail(k);
                if constexpr (Unit)
                    corner += (-trig.col(k).head(k)) * head.transpose();
                else
                    corner += divide(-trig.col(k).head(k), trig.diag().head(k)) * head.transpose();
            }
        }
        else {
            for (size_t i = 0; i + 1 < getOrder(); ++i) {
                auto corner = target.bottomLeftCorner(i + 1, i + 1);
                auto head = target.row(i).head(i + 1);
                if constexpr (Unit)
                    corner += (-trig.col(i).tail(i + 1)) * head.transpose();
                else
                    corner += divide(-trig.col(i).tail(i + 1), trig.diag().tail(i + 1)) * head.transpose();
            }
        }
    }

    template<Matrix M> requires(instanceof_tx<M, MatrixTrig>)
    auto&& Inverse<M>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.trig);
    }

    template<Matrix M> requires(instanceof_tx<M, MatrixTrig>)
    __host__ __device__ consteval int Inverse<M>::getMajor() noexcept {
        return std::remove_cvref_t<M>::getMajor();
    }
}

#ifdef PHYSICA_MKL
    #include "Inverse_MKL.h"
#endif
