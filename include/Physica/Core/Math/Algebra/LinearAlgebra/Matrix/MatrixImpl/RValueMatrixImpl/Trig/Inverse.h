/*
 * Copyright 2026 Weibo He.
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
    template<Matrix M> requires(instanceof_tx<MatrixTrig, M>)
    class Inverse<M> : public RValueMatrix<Inverse<M>> {
        using This = Inverse<M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
        using typename Base::Trv;
        using typename Base::Tc;
        using typename Base::Tm;
    private:
        LazyDestroy<M&&> trig;
    public:
        explicit Inverse(M trig);
        Inverse(const This&) = default;
        Inverse(This&&) noexcept = default;
        ~Inverse() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Operations */
        void assign(Matrix auto& target) const;
        void assign_base(Matrix auto& target) const noexcept;
        void assign_mkl(Matrix auto& target) const;
        /* Getters */
        [[nodiscard]] const auto& getExpr() const noexcept { return trig; }
        [[nodiscard]] size_t getRow() const noexcept { return trig.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return getRow(); }
    };

    template<Matrix M> requires(instanceof_tx<MatrixTrig, M>)
    Inverse<M>::Inverse(M trig) : trig(std::forward<M>(trig)) {}

    template<Matrix M> requires(instanceof_tx<MatrixTrig, M>)
    void Inverse<M>::assign(Matrix auto& target) const {
        using Expr = Traits<M>::ExprType;
        if constexpr (Internal::EnableMKL<Expr, decltype(target)>::value)
            assign_mkl(target);
        else
            assign_base(target);
    }

    template<Matrix M> requires(instanceof_tx<MatrixTrig, M>)
    void Inverse<M>::assign_base(Matrix auto& target) const noexcept {
        target.assert_assign(*this);
        if constexpr (!Traits<M>::Upper && Traits<M>::Unit) {
            (-trig).assign(target);
            target.diag() = Trv(1);
            for (size_t i = 1; i < getCol() - 1; ++i) {
                auto corner = target.bottomLeftCorner(i + 1, i);

                auto row = target.row(i);
                auto col = target.col(i);
                auto head = row.head(i);
                auto tail = col.tail(i + 1);
                corner += tail * head.transpose();
            }
        }
        else
            noImpl(__func__);
    }
}

#ifdef PHYSICA_MKL
    #include "Inverse_MKL.h"
#endif
