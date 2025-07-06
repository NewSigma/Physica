/*
 * Copyright 2024-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/UnitVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RValueMatrix.h"

namespace Physica {
    template<Matrix M, Vector V> class MatExpVecProd;

    template<Matrix M>
    class MatrixExp : public RValueMatrix<MatrixExp<M>> {
        using This = MatrixExp<M>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
        using typename Base::Tr;
    private:
        const LazyDestroy<M> m;
    public:
        MatrixExp(M&& m_);
        MatrixExp(const This&) = default;
        MatrixExp(This&&) noexcept = default;
        ~MatrixExp() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;

        template<Vector V>
        [[nodiscard]] auto operator*(V&& v) const& noexcept;
        template<Vector V>
        [[nodiscard]] auto operator*(V&& v) && noexcept;
        /* Operations */
        void assign(Matrix auto& target) const;
        [[nodiscard]] T calc(size_t, size_t) const { noImpl("calc() is low performance and should be avoided"); }

        [[nodiscard]] Tr calcTraceMu() const;
        /* Getters */
        [[nodiscard]] const auto& getMatrix() const noexcept { return m; }
        [[nodiscard]] size_t getRow() const noexcept { return m.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return m.getCol(); }
    };

    template<Matrix M>
    MatrixExp<M>::MatrixExp(M&& m_) : m(std::forward<M>(m_)) {
        assert(m.getRow() == m.getCol());
    }

    template<Matrix M>
    template<Vector V>
    auto MatrixExp<M>::operator*(V&& v) const& noexcept {
        return MatExpVecProd<const This&, V&&>(*this, std::forward<V>(v));
    }

    template<Matrix M>
    template<Vector V>
    auto MatrixExp<M>::operator*(V&& v) && noexcept {
        return MatExpVecProd<This&&, V&&>(std::move(*this), std::forward<V>(v));
    }

    template<Matrix M>
    void MatrixExp<M>::assign(Matrix auto& target) const {
        const Tr traceMu = calcTraceMu();
        const auto params = ((*this) * VectorND<T>(getRow())).template calcParam<Sequential>(traceMu);
        for (size_t i = 0; i < getCol(); ++i) {
            auto col = target.col(i);
            ((*this) * UnitVector<T>(i, getRow())).assign(col, traceMu, params);
        }
    }

    template<Matrix M>
    auto MatrixExp<M>::calcTraceMu() const -> Tr {
        const T trace = m.trace();
        assert(trace.imag().isZero() && "[Error]: Not implemented");
        return trace.real() / Tr(getRow());
    }

    template<Matrix M>
    [[nodiscard]] inline auto exp(M&& m) noexcept {
        return MatrixExp<M&&>(std::forward<M>(m));
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<MatrixExp<M>> : public Traits<std::remove_cvref_t<M>> {
    public:
        constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::AnyStorage;
    };
}

#include "MatExpVecProd.h"
