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

#include "Physica/Core/Scalar/Scalar.h"
#include "LatticeModel.h"

namespace Physica {
    template<Scalar T, int Dim>
    class TransIsing : public LatticeModel<Dim> {
        static_assert(Dim == 1, "Not implemented");
        using This = TransIsing<T, Dim>;
        using Base = LatticeModel<Dim>;
    public:
        using ScalarType = T;
    private:
        T couplingJ;
        T transG;
    public:
        TransIsing(Base base, T couplingJ_, T transG_);
        TransIsing(const This&) = default;
        TransIsing(This&&) noexcept = default;
        ~TransIsing() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] T getCouplingJ() const noexcept { return couplingJ; }
        [[nodiscard]] T getTransG() const noexcept { return transG; }
    };

    template<Scalar T, int Dim>
    TransIsing<T, Dim>::TransIsing(Base base, T couplingJ_, T transG_) : Base(std::move(base)), couplingJ(couplingJ_), transG(transG_) {}

    template<Scalar T, int Dim>
    void TransIsing<T, Dim>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        couplingJ.swap(obj.couplingJ);
        transG.swap(obj.transG);
    }
}

namespace Physica {
    template<Scalar T, int Dim_>
    class Traits<TransIsing<T, Dim_>> {
    public:
        constexpr static int Dim = Dim_;
        constexpr static int SiteDOF = 2;
    };
}
