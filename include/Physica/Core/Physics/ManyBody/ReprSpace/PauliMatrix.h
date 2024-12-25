/*
 * Copyright 2024 Weibo He.
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

namespace Physica::Core {
    enum class PauliIndex {
        I,
        X,
        Y,
        Z
    };

    template<Scalar T, PauliIndex Idx>
    class PauliMatrix : public RValueMatrix<PauliMatrix<T, Idx>> {
        using This = PauliMatrix<T, Idx>;
        using ElemType = std::conditional<Idx == PauliIndex::Y, typename T::ComplexType, T>::type;

        int site;
    public:
        PauliMatrix(int site_);
        PauliMatrix(const This&) = default;
        PauliMatrix(This&&) noexcept = default;
        ~PauliMatrix() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        constexpr static ElemType calc(size_t row, size_t col) noexcept;

        template<int Dim, int NumSite>
        ElemType apply(SpinState<Dim, NumSite>& psi) const noexcept;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        constexpr static size_t getRow() noexcept { return 2; }
        constexpr static size_t getCol() noexcept { return 2; }
    };

    template<Scalar T, PauliIndex Idx>
    PauliMatrix<T, Idx>::PauliMatrix(int site_) : site(site_) {}

    template<Scalar T, PauliIndex Idx>
    constexpr auto PauliMatrix<T, Idx>::calc(size_t row, size_t col) noexcept -> ElemType {
        using enum PauliIndex;
        if constexpr (Idx == I)
            return ElemType(row == col ? 1 : 0);
        else if constexpr (Idx == X)
            return ElemType(row == col ? 0 : 1);
        else if constexpr (Idx == Y)
            return ElemType(0, row == col ? 0 : (row == 0 ? -1 : 1));
        else
            return ElemType(row == col ? (row == 0 ? 1 : -1) : 0);
    }

    template<Scalar T, PauliIndex Idx>
    template<int Dim, int NumSite>
    auto PauliMatrix<T, Idx>::apply(SpinState<Dim, NumSite>& psi) const noexcept -> ElemType {
        using enum PauliIndex;
        if constexpr (Idx == I)
            return ElemType(1);
        else if constexpr (Idx == X) {
            psi = psi.flip(site);
            return ElemType(1);
        }
        else if constexpr (Idx == Y) {
            const bool flag = psi.isSpinUp(site);
            psi = psi.flip(site);
            return ElemType(0, flag ? 1 : -1);
        }
        else {
            return ElemType(psi.isSpinUp(site) ? 1 : -1);
        }
    }

    template<Scalar T, PauliIndex Idx>
    void PauliMatrix<T, Idx>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(site, obj.site);
    }
}

namespace Physica {
    template<Scalar T, PauliIndex I>
    class Traits<PauliMatrix<T, I>> {
        static_assert(!T::isComplex);
    public:
        using ScalarType = T;
        constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::AnyStorage;
        constexpr static size_t RowAtCompile = 2;
        constexpr static size_t ColAtCompile = 2;
        constexpr static size_t SizeAtCompile = 4;
    };
}
