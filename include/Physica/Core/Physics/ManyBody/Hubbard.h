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

#include <Physica/Core/Math/Transform/FFT.h>
#include "LatticeHamilton.h"
#include "ReprSpace/SpinRepr.h"

namespace Physica::Core {
    template<class ScalarType, class ReprType> class Hubbard;

    namespace Internal {
        template<class T, class U>
        class Traits<Hubbard<T, U>>
                : public Traits<LatticeHamilton<Hubbard<T, U>>> {
        public:
            using ScalarType = T;
            using ReprType = U;
        };
    }
    /**
     * Refer to [1] for applied symmetries
     * 
     * Reference:
     * [1] J. Korean Phys. Soc. 76, 670–683 (2020); https://doi.org/10.3938/jkps.76.670
     */
    template<class ScalarType, class ReprType>
    class Hubbard : public LatticeHamilton<Hubbard<ScalarType, ReprType>> {
        using This = Hubbard<ScalarType, ReprType>;
        using Base = LatticeHamilton<This>;
        using RealType = typename ScalarType::RealType;
        using VectorType = Vector<ScalarType>;
        using FFTType = FFT<RealType, 1>;
        using typename Base::StateType;
        using typename Base::IndexType;

        constexpr static bool IsTransInvariant = Internal::Traits<ReprType>::IsTransInvariant;
        static_assert(std::is_base_of<ReprBasis<ReprType>, ReprType>::value, "[Error]: ReprType is not a representation");
        static_assert((IsTransInvariant && Base::isComplex) || !IsTransInvariant, "[Error]: Use complex scalar if translational invariance is enabled");
        static_assert(!IsTransInvariant || (Base::Dim == 1), "[Error]: Trans invariantce is not implemented in high dimension");
    public:
        using typename Base::DimArray;
        using Base::Dim;
        using Base::NumSite;
    private:
        ReprType repr;
        RealType hoppingT;
        RealType repelU;
        FFTType planProvider;
    public:
        Hubbard() = default;
        Hubbard(DimArray superSize, unsigned int numUnitCellSite, ReprType repr_, RealType hoppingT_, RealType repelU_);
        Hubbard(const This&) = default;
        Hubbard(This&&) noexcept = default;
        ~Hubbard() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        template<class AnyVector>
        [[nodiscard]] Vector<ScalarType> operator*(const RValueVector<AnyVector>& v) const;
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const ReprType& getRepr() const noexcept { return repr; }
        [[nodiscard]] size_t getNumState() const noexcept { return repr.getNumState(); }
        [[nodiscard]] RealType getHoppingT() const noexcept { return hoppingT; }
        [[nodiscard]] RealType getRepelU() const noexcept { return repelU; }
    protected:
        RealType repelElem(StateType psi) const;
        RealType hoppingElem(StateType rowPsi, StateType colPsi) const;
        void sumHopping(VectorType& target, ScalarType value, StateType psi) const noexcept;
        void sumHopping(VectorType& target, FFTType& fft, ScalarType factor, StateType psi) const;
    private:
        void dotImpl1D(Vector<ScalarType>& result, ScalarType factor, size_t index) const;
        void dotImplND(Vector<ScalarType>& result, ScalarType factor, size_t index) const;
    };
}

#include "HubbardImpl.h"
