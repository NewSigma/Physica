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
#include "ReprSpace/SpinRepr.h"

namespace Physica::Core {
    template<class ScalarType, class ReprType> class Hubbard1D;

    namespace Internal {
        template<class T, class U>
        class Traits<Hubbard1D<T, U>>
                : public Traits<LatticeHamilton<Hubbard1D<T, U>>> {
        public:
            using ScalarType = T;
            using ReprType = U;
        };
    }
    /**
     * Refer to [1] for applied symmetries
     * 
     * Reference:
     * [1] Computers in Physics 7, 400 (1993); https://doi.org/10.1063/1.4823192
     */
    template<class ScalarType, class ReprType>
    class Hubbard1D : public LatticeHamilton<Hubbard1D<ScalarType, ReprType>> {
        static_assert(std::is_base_of<ReprBasis<ReprType>, ReprType>::value, "[Error]: Invalid representation");
        using This = Hubbard1D<ScalarType, ReprType>;
        using Base = LatticeHamilton<This>;
        using typename Base::LatticeType;
        using typename Base::StateType;
    private:
        ReprType repr;
        ScalarType hoppingT;
        ScalarType repelU;
    public:
        Hubbard1D() = default;
        Hubbard1D(LatticeType lattice, ReprType repr_, ScalarType hoppingT_, ScalarType repelU_);
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
        [[nodiscard]] const ReprType& getRepr() const noexcept { return repr; }
        [[nodiscard]] size_t getNumState() const noexcept { return repr.getNumState(); }
        [[nodiscard]] ScalarType getHoppingT() const noexcept { return hoppingT; }
        [[nodiscard]] ScalarType getRepelU() const noexcept { return repelU; }
        [[nodiscard]] size_t getNumSuperCellSite() const noexcept { return Base::getLattice().getNumSuperCellSite(); }
    };
}

#include "Hubbard1DImpl.h"
