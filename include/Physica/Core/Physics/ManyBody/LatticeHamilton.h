/*
 * Copyright 2024 WeiBo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "PeriodicLattice.h"

namespace Physica::Core {
    template<class Derived> class LatticeHamilton;

    namespace Internal {
        template<class Derived>
        class Traits<LatticeHamilton<Derived>> {
        public:
            constexpr static int Option = MatrixOption::AnyMajor | MatrixOption::Element;
            constexpr static size_t RowAtCompile = Utils::Dynamic;
            constexpr static size_t ColumnAtCompile = Utils::Dynamic;
            constexpr static size_t MaxRowAtCompile = Utils::Dynamic;
            constexpr static size_t MaxColumnAtCompile = Utils::Dynamic;
            constexpr static size_t SizeAtCompile = Utils::Dynamic;
            constexpr static size_t MaxSizeAtCompile = Utils::Dynamic;
        };
    }

    template<class Derived>
    class LatticeHamilton : public RValueMatrix<Derived> {
        using This = LatticeHamilton<Derived>;
        using Base = RValueMatrix<Derived>;
        using DerivedTraits = Internal::Traits<Derived>;
        constexpr static unsigned int Dim = DerivedTraits::Dim;
    public:
        using ScalarType = typename DerivedTraits::ScalarType;
        using StateType = typename DerivedTraits::StateType;
        using LatticeType = PeriodicLattice<Dim>;
    private:
        LatticeType lattice;
    public:
        LatticeHamilton(LatticeType lattice_);
        LatticeHamilton(const LatticeHamilton&) = default;
        LatticeHamilton(LatticeHamilton&&) noexcept = default;
        ~LatticeHamilton() = default;
        /* Operators */
        LatticeHamilton& operator=(LatticeHamilton obj) noexcept { swap(obj); return *this; }
        template<class VectorType>
        [[nodiscard]] Vector<ScalarType> operator*(const RValueVector<VectorType>& v) const;
        /* Operations */
        [[nodiscard]] size_t stateToIndex(StateType state) const noexcept { return Base::getDerived().stateToIndex(state); }
        [[nodiscard]] StateType indexToState(size_t index) const noexcept { return Base::getDerived().indexToState(index); }
        [[nodiscard]] const This& hermite() const noexcept { return *this; }
        void swap(LatticeHamilton& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const LatticeType& getLattice() const noexcept { return lattice; }
        [[nodiscard]] size_t getNumState() const noexcept { return Base::getDerived().getNumState(); }
        [[nodiscard]] size_t getRow() const noexcept { return getNumState(); }
        [[nodiscard]] size_t getColumn() const noexcept { return getNumState(); }
    };

    template<class Derived>
    LatticeHamilton<Derived>::LatticeHamilton(LatticeType lattice_) : lattice(std::move(lattice_)) {}

    template<class Derived>
    template<class VectorType>
    Vector<typename LatticeHamilton<Derived>::ScalarType> LatticeHamilton<Derived>::operator*(const RValueVector<VectorType>& v) const {
        assert(getColumn() == v.getLength() && "[Error]: Dimensions do not match");
        return Base::getDerived() * v;
    }

    template<class Derived>
    void LatticeHamilton<Derived>::swap(LatticeHamilton& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        lattice.swap(obj.lattice);
    }
}
