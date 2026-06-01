/*
 * Copyright 2025-2026 Weibo He.
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

#include "MatrixTrig.h"

namespace Physica {
    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig>)
    class GEMM<M1, M2> : public RValueMatrix<GEMM<M1, M2>> {
        using This = GEMM<M1, M2>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
        using typename Base::Tc;
    private:
        LazyDestroy<M1> trig;
        LazyDestroy<M2> rhs;
    public:
        GEMM(M1&& trig, M2&& rhs);
        GEMM(const This&) = default;
        GEMM(This&&) noexcept = default;
        ~GEMM() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto& target) const;
        void assign_mkl(Matrix auto& target) const noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const { return trig.getRow(); }
        [[nodiscard]] size_t getCol() const { return rhs.getCol(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
    };

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig>)
    GEMM<M1, M2>::GEMM(M1&& trig, M2&& rhs) : trig(std::forward<M1>(trig)), rhs(std::forward<M2>(rhs)) {}

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig>)
    void GEMM<M1, M2>::assign(Matrix auto& target) const {
        if constexpr (HasMKL())
            assign_mkl(target);
        else
            noImpl();
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig>)
    auto&& GEMM<M1, M2>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M1>(self.trig);
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig>)
    auto&& GEMM<M1, M2>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M2>(self.rhs);
    }
}

#ifdef PHYSICA_MKL
    #include "GEMM_MKL.h"
#endif
