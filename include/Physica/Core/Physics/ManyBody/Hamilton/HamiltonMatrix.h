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
#include "Physica/Core/Physics/ManyBody/Model/LatticeModel.h"

namespace Physica {
    template<class Derived>
    class HamiltonMatrix : public RValueMatrix<Derived> {
        using This = HamiltonMatrix<Derived>;
        using Base = RValueMatrix<Derived>;

        using T = Traits<Derived>::ScalarType;
        using ReprType = Traits<Derived>::ReprType;
        constexpr static BoundaryCond Boundary = Traits<Derived>::Boundary;
    protected:
        using StateType = ReprType::StateType;
    public:
        constexpr static unsigned int Dim = ReprType::Dim;
        constexpr static unsigned int NumSite = StateType::NumSite;
        constexpr static unsigned int SiteDOF = StateType::SiteDOF;
        constexpr static bool IsTransInvariant = Traits<ReprType>::IsTransInvariant;

        static_assert((IsTransInvariant && T::isComplex()) || !IsTransInvariant, "[Error]: Use complex scalar if translational invariance is enabled");
        static_assert(!IsTransInvariant || (Dim == 1), "[Error]: Translational invariance is not implemented in high dimension");
        static_assert((Boundary == BoundaryCond::PBC) || !IsTransInvariant, "[Error]: Translational invariance is not implemented for other boundary conditions");
        static_assert((Boundary != BoundaryCond::OBC) || !IsTransInvariant, "[Error]: Translational invariance does not support OBC");
        static_assert((Boundary != BoundaryCond::TBC) || T::isComplex(), "[Error]: Twisted boundary condition gives complex hamilton");
    public:
        ~HamiltonMatrix() = default;
        /* Operations */
        [[nodiscard]] decltype(auto) transpose(this auto&&) noexcept;
        [[nodiscard]] const Derived& hermite() const noexcept { return Base::getDerived(); }
        /* Getters */
        [[nodiscard]] const auto& getModel() const noexcept { return Base::getDerived().getModel(); }
        [[nodiscard]] const auto& getRepr() const noexcept { return Base::getDerived().getRepr(); }
        [[nodiscard]] size_t getNumState() const noexcept { return getRepr().getNumState(); }
        [[nodiscard]] size_t getRow() const noexcept { return getNumState(); }
        [[nodiscard]] size_t getCol() const noexcept { return getNumState(); }
    protected:
        HamiltonMatrix() = default;
        HamiltonMatrix(const This&) = default;
        HamiltonMatrix(This&&) noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
    };

    template<class Derived>
    decltype(auto) HamiltonMatrix<Derived>::transpose(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (T::isComplex()) {
            using X = Base; // FIXME: clang 22 rejects valid
            [[maybe_unused]] auto x = sizeof(X);
            return std::forward<Self>(self).X::transpose();
        }
        else
            return std::forward<Self>(self);
    }
}

namespace Physica {
    template<class Derived>
    class Traits<HamiltonMatrix<Derived>> {
    public:
        constexpr static int Major = MatrixMajor::BothMajor;
    };
}
