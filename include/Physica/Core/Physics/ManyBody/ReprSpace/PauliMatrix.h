/*
 * Copyright 2024-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RValueMatrix.h"
#include "State/SpinState.h"

namespace Physica {
    enum class PauliIndex : char {
        I = 0,
        X = 1,
        Y = 2,
        Z = 3
    };

    template<Scalar U, PauliIndex Idx>
    class PauliMatrix : public RValueMatrix<PauliMatrix<U, Idx>> {
        using This = PauliMatrix<U, Idx>;
        using Base = RValueMatrix<This>;
    protected:
        using typename Base::T;
    private:
        int site = 0;
    public:
        PauliMatrix() = default;
        PauliMatrix(int site_);
        PauliMatrix(const This&) = default;
        PauliMatrix(This&&) noexcept = default;
        ~PauliMatrix() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        constexpr static T calc(size_t row, size_t col) noexcept;

        template<int Dim, int NumSite>
        T apply(SpinState<Dim, NumSite>& psi) const noexcept;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        constexpr static size_t getRow() noexcept { return 2; }
        constexpr static size_t getCol() noexcept { return 2; }
    };

    template<Scalar U, PauliIndex Idx>
    PauliMatrix<U, Idx>::PauliMatrix(int site_) : site(site_) {}

    template<Scalar U, PauliIndex Idx>
    constexpr auto PauliMatrix<U, Idx>::calc(size_t row, size_t col) noexcept -> T {
        using enum PauliIndex;
        if constexpr (Idx == I)
            return T(row == col ? 1 : 0);
        else if constexpr (Idx == X)
            return T(row == col ? 0 : 1);
        else if constexpr (Idx == Y)
            return T(0, row == col ? 0 : (row == 0 ? -1 : 1));
        else
            return T(row == col ? (row == 0 ? 1 : -1) : 0);
    }

    template<Scalar U, PauliIndex Idx>
    template<int Dim, int NumSite>
    auto PauliMatrix<U, Idx>::apply(SpinState<Dim, NumSite>& psi) const noexcept -> T {
        using enum PauliIndex;
        if constexpr (Idx == I)
            return T(1);
        else if constexpr (Idx == X) {
            psi = psi.flip(site);
            return T(1);
        }
        else if constexpr (Idx == Y) {
            const bool flag = psi.isSpinUp(site);
            psi = psi.flip(site);
            return T(0, flag ? 1 : -1);
        }
        else {
            return T(psi.isSpinUp(site) ? 1 : -1);
        }
    }

    template<Scalar U, PauliIndex Idx>
    void PauliMatrix<U, Idx>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(site, obj.site);
    }
}

namespace Physica {
    template<Scalar U, PauliIndex Idx>
    class Traits<PauliMatrix<U, Idx>> {
        static_assert(!U::isComplex);
    public:
        using ScalarType = std::conditional<Idx == PauliIndex::Y, typename U::ComplexType, U>::type;
        constexpr static int Option = MatrixMajor::BothMajor;
        constexpr static size_t RowAtCompile = 2;
        constexpr static size_t ColAtCompile = 2;
        constexpr static size_t SizeAtCompile = 4;
    };
}
