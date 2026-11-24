/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
        [[nodiscard]] size_t getOrder() const noexcept { return getNumState(); }
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept { return true; }
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return MatrixMajor::BothMajor; }
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
