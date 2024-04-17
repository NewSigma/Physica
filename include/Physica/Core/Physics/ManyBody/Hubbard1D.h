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

#include "LatticeHamilton.h"
#include "State/SpinElectron.h"

namespace Physica::Core {
    template<class ScalarType> class Hubbard1D;

    namespace Internal {
        template<class T>
        class Traits<Hubbard1D<T>> : public Traits<LatticeHamilton<Hubbard1D<T>>> {
        public:
            using ScalarType = T;
            constexpr static unsigned int Dim = 1;
            constexpr static unsigned int SiteDOF = 4;
        };
    }
    /**
     * Reference:
     * [1] Computers in Physics 7, 400 (1993); https://doi.org/10.1063/1.4823192
     */
    template<class ScalarType>
    class Hubbard1D : public LatticeHamilton<Hubbard1D<ScalarType>> {
        using This = Hubbard1D<ScalarType>;
        using Base = LatticeHamilton<This>;
        using LatticeType = typename Base::LatticeType;
        using StateType = SpinElectron;
    private:
        ScalarType hoppingT;
        ScalarType repelU;
        unsigned int numSpinUp;
        unsigned int numSpinDown;
    public:
        Hubbard1D(LatticeType lattice, ScalarType hoppingT_, ScalarType repelU_, unsigned int numSpinUp_, unsigned int numSpinDown_);
        Hubbard1D(const Hubbard1D&) = default;
        Hubbard1D(Hubbard1D&&) noexcept = default;
        ~Hubbard1D() = default;
        /* Operators */
        Hubbard1D& operator=(Hubbard1D obj) noexcept { swap(obj); return *this; }
        template<class VectorType>
        [[nodiscard]] Vector<ScalarType> operator*(const RValueVector<VectorType>& v) const;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const;
        void swap(Hubbard1D& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumSuperCellSite() const noexcept { return Base::getLattice().getNumSuperCellSite(); }
    private:
        [[nodiscard]] size_t stateToIndex(StateType state) const noexcept;
        [[nodiscard]] StateType indexToState(size_t index) const noexcept;
        void stateAdd(Vector<ScalarType>& target, StateType state, ScalarType value) const noexcept;
        void checkState(StateType state) const noexcept;
    };
}

#include "Hubbard1DImpl.h"
