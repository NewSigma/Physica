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

#include "Physica/Core/Physics/ManyBody/ReprSpace/ReprBasis.h"
#include "Physica/Core/Physics/ManyBody/Model/TransIsing.h"
#include "HamiltonMatrix.h"

namespace Physica::Core {
    template<Scalar T, Representation U>
    class TransIsingMatrix : public HamiltonMatrix<TransIsingMatrix<T, U>>, public TransIsing<T, U::Dim> {
        using This = TransIsingMatrix<T, U>;
        using Base = HamiltonMatrix<This>;
        using ModelBase = TransIsing<T, U::Dim>;

        using StateType = U::StateType;
    public:
        constexpr static int Dim = U::Dim;
        constexpr static int NumSite = StateType::NumSite;
        constexpr static int SiteDOF = StateType::SiteDOF;
    private:
        U repr;
    public:
        TransIsingMatrix() = default;
        TransIsingMatrix(ModelBase model, U repr_);
        TransIsingMatrix(const This&) = default;
        TransIsingMatrix(This&&) noexcept = default;
        ~TransIsingMatrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const ModelBase& getModel() const noexcept { return *this; }
        [[nodiscard]] const U& getRepr() const noexcept { return repr; }
        using ModelBase::getCouplingJ;
        using ModelBase::getTransG;
    };

    template<Scalar T, Representation U>
    TransIsingMatrix<T, U>::TransIsingMatrix(ModelBase model, U repr_) : ModelBase(std::move(model)), repr(std::move(repr_)) {}

    template<Scalar T, Representation U>
    T TransIsingMatrix<T, U>::calc(size_t row, size_t col) const {
        const auto psi1 = repr[row];
        const auto mask = StateType::makeFullMask() >> 1;
        if (row == col) {
            const auto bits = psi1.getOccupyBits();
            const int numAntiSpin = std::popcount((bits >> 1) ^ (bits & mask));
            assert(numAntiSpin < NumSite);
            return -getCouplingJ() * T((int(NumSite - 1) - numAntiSpin * 2));
        }

        const auto psi2 = repr[col];
        if ((psi1 ^ psi2).getNumParticle() == 1)
            return -getTransG();
        return 0;
    }

    template<Scalar T, Representation U>
    void TransIsingMatrix<T, U>::swap(This& __restrict obj) noexcept {
        ModelBase::swap(obj);
        repr.swap(obj.repr);
    }
}

namespace Physica {
    template<Scalar T, Representation U>
    class Traits<TransIsingMatrix<T, U>> : public Traits<HamiltonMatrix<TransIsingMatrix<T, U>>> {
    public:
        using ScalarType = T;
        using ReprType = U;
    };
}

#include "TransIsingVecProduct.h" // IWYU pragma: export
