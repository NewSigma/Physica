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

#include "EmptyForceModel.h"

namespace Physica {
    /**
     * Reference:
     * [1] Phys. Rev. Lett. 134, 246401; https://doi.org/10.1103/9gkl-w2lm
     */
    template<Scalar T>
    class DoubleWell : private EmptyForceModel<T, 1> {
        using Base = EmptyForceModel<T, 1>;
        using This = DoubleWell<T>;
    public:
        using typename Base::MDCellType;
        using typename Base::LatticeMatrix;
        using typename Base::ForceConstMatrix;
        using PositionMatrix = MDCellType::PositionMatrix;
    private:
        T strength;
        T widthR;
    public:
        DoubleWell(T strength_, T widthR_);
        DoubleWell(const This&) = default;
        DoubleWell(This&&) noexcept = default;
        ~DoubleWell() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] T potentialV(const MDCellType& cell) const;

        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force(const MDCellType& cell) const;
        template<ExecutePolicy P>
        void forceAsync(const MDCellType& cell, Vector auto& result) const;
        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force_short(const MDCellType& cell) const { return force<P>(cell); }
        using Base::force_long;

        [[nodiscard]] T forceConst([[maybe_unused]] const MDCellType& cell, size_t dof1, size_t dof2) const;
        [[nodiscard]] ForceConstMatrix forceConst(const MDCellType& cell) const;

        [[nodiscard]] LatticeMatrix virial(const MDCellType& cell) const;
        void swap(This& __restrict obj) noexcept;
        /* Static members */
        [[nodiscard]] consteval static bool isPeriodBoundary() noexcept { return false; }
    };

    template<Scalar T>
    DoubleWell<T>::DoubleWell(T strength_, T widthR_)
            : strength(strength_)
            , widthR(widthR_) {
        assert(strength.isPositive() && widthR.isPositive());
    }

    template<Scalar T>
    T DoubleWell<T>::potentialV(const MDCellType& cell) const {
        return square_elem(square_elem(cell.getPos()) - square(widthR)).sum() * strength;
    }

    template<Scalar T>
    template<ExecutePolicy P>
    VectorND<T> DoubleWell<T>::force(const MDCellType& cell) const {
        VectorND<T> result(cell.getDOF());
        forceAsync<P>(cell, result);
        return result;
    }

    template<Scalar T>
    template<ExecutePolicy P>
    void DoubleWell<T>::forceAsync(const MDCellType& cell, Vector auto& result) const {
        auto x = cell.getPos().col(0);
        result = T(-4.0) * strength * hadamard(x, square(x) - square(widthR));
    }

    template<Scalar T>
    T DoubleWell<T>::forceConst([[maybe_unused]] const MDCellType& cell, size_t dof1, size_t dof2) const {
        noImpl();
    }

    template<Scalar T>
    auto DoubleWell<T>::forceConst(const MDCellType& cell) const -> ForceConstMatrix {
        noImpl();
    }

    template<Scalar T>
    auto DoubleWell<T>::virial(const MDCellType& cell) const -> LatticeMatrix {
        noImpl();
    }

    template<Scalar T>
    void DoubleWell<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        strength.swap(obj.strength);
        widthR.swap(obj.widthR);
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<DoubleWell<T>> {
    public:
        constexpr static bool IsContractable = false;
    };
}
