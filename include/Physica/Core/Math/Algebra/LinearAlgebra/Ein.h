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

#include <cassert>
#include <utility>
#include "Physica/Core/Math/Algebra/LinearAlgebra/IndexVar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/Tensor.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Core/Utils/Container/Array.h"
#include "Physica/Core/Utils/MetaProgramming.h"

namespace Physica {
    template<class LHS, class RHS> class Contract;

    template<class Expr>
    class Ein {
        using This = Ein<Expr>;
    public:
        constexpr static int NDim = []() consteval static noexcept {
            if constexpr (Tensor<Expr>)
                return std::remove_cvref_t<Expr>::NDim;
            else if constexpr (Matrix<Expr>)
                return 2;
            else {
                static_assert(Vector<Expr>, "[Error]: Einstein notation requires a tensor, matrix or vector");
                return 1;
            }
        }();
        using IndexType = Array<size_t, NDim>;
    protected:
        using T = std::remove_cvref_t<Expr>::ScalarType;
    public:
        using ScalarType = T;
    private:
        decay_rvalue_t<Expr> expr;
        Array<const Var*, NDim> vars;
    public:
        Ein(Expr&& expr_, Array<const Var*, NDim> vars_) noexcept;
        Ein(const This&) = default;
        Ein(This&&) noexcept = default;
        ~Ein() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        This& operator=(const instanceof<Contract> auto& expr);

        [[nodiscard]] decltype(auto) operator[](const IndexType& indices);
        [[nodiscard]] decltype(auto) operator[](std::same_as<size_t> auto... indices);
        [[nodiscard]] auto operator*(this auto&& lhs, instanceof<Ein> auto&& rhs) noexcept;
        /* Operations */
        [[nodiscard]] T calc(const IndexType& indices) const;
        [[nodiscard]] T calc(std::same_as<size_t> auto... indices) const;
        /* Getters */
        [[nodiscard]] auto&& getExpr(this auto&&) noexcept;
        [[nodiscard]] const Array<const Var*, NDim>& getVars() const noexcept { return vars; }
        [[nodiscard]] IndexType getShape() const noexcept;
        [[nodiscard]] size_t dim(int index) const noexcept;
    };

    template<class Expr>
    Ein<Expr>::Ein(Expr&& expr_, Array<const Var*, NDim> vars_) noexcept
            : expr(std::forward<Expr>(expr_)), vars(std::move(vars_)) {}

    template<class Expr>
    auto Ein<Expr>::operator=(const instanceof<Contract> auto& expr) -> This& {
        expr.assign(*this);
        return *this;
    }

    template<class Expr>
    decltype(auto) Ein<Expr>::operator[](const IndexType& indices) {
        if constexpr (Vector<Expr>)
            return getExpr()[indices[0]];
        else if constexpr (Matrix<Expr>)
            return getExpr()[indices[0], indices[1]];
        else
            return getExpr()[indices];
    }

    template<class Expr>
    decltype(auto) Ein<Expr>::operator[](std::same_as<size_t> auto... indices) {
        return getExpr()[indices...];
    }

    template<class Expr>
    auto Ein<Expr>::operator*(this auto&& lhs, instanceof<Ein> auto&& rhs) noexcept {
        return Contract(std::forward<decltype(lhs)>(lhs), std::forward<decltype(rhs)>(rhs));
    }

    template<class Expr>
    auto Ein<Expr>::calc(const IndexType& indices) const -> T {
        if constexpr (Vector<Expr>)
            return getExpr().calc(indices[0]);
        else if constexpr (Matrix<Expr>)
            return getExpr().calc(indices[0], indices[1]);
        else
            return getExpr().calc(indices);
    }

    template<class Expr>
    auto Ein<Expr>::calc(std::same_as<size_t> auto... indices) const -> T {
        return getExpr().calc(indices...);
    }

    template<class Expr>
    auto&& Ein<Expr>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Expr>(self.expr);
    }

    template<class Expr>
    auto Ein<Expr>::getShape() const noexcept -> IndexType {
        if constexpr (Vector<Expr>)
            return IndexType{expr.getLength()};
        else if constexpr (Matrix<Expr>)
            return IndexType{expr.getRow(), expr.getCol()};
        else
            return expr.getShape();
    }

    template<class Expr>
    size_t Ein<Expr>::dim(int index) const noexcept {
        if constexpr (Vector<Expr>) {
            assert(index == 0 && "[Error]: A vector has exactly 1 dimension");
            return expr.getLength();
        }
        else if constexpr (Matrix<Expr>) {
            assert((index == 0 || index == 1) && "[Error]: A matrix has exactly 2 dimensions");
            return index == 0 ? expr.getRow() : expr.getCol();
        }
        else
            return expr.dim(index);
    }
}

#include "Contract.h"
