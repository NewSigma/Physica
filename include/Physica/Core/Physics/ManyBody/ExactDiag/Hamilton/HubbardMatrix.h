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
#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/SparseVector.h>
#include <Physica/Core/Physics/ManyBody/ExactDiag/ReprSpace/ReprImpl/ReprSpace.h>
#include <Physica/Core/Physics/ManyBody/Model/Hubbard.h>
#include "HamiltonMatrix.h"

namespace Physica::Core {
    class ThreadExecutor;
    /**
     * Refer to [1] for applied symmetries
     * 
     * Reference:
     * [1] J. Korean Phys. Soc. 76, 670–683 (2020); https://doi.org/10.3938/jkps.76.670
     */
    template<class ScalarType, class ReprType>
    class HubbardMatrix
            : public HamiltonMatrix<HubbardMatrix<ScalarType, ReprType>>
            , public Hubbard<typename ScalarType::RealType, ReprType::Dim> {
        using RealType = typename ScalarType::RealType;
        using This = HubbardMatrix<ScalarType, ReprType>;
        using Base = HamiltonMatrix<This>;
        using ModelBase = Hubbard<RealType, ReprType::Dim>;
        using FFTType = FFT<RealType, 1>;
        using StateType = typename ReprType::StateType;
        using typename ModelBase::IndexType;
        constexpr static unsigned int SiteDOF = StateType::SiteDOF;
        constexpr static bool IsTransInvariant = Traits<ReprType>::IsTransInvariant;
    public:
        constexpr static unsigned int Dim = ReprType::Dim;
        constexpr static unsigned int NumSite = StateType::NumSite;
        static_assert((IsTransInvariant && ScalarType::isComplex) || !IsTransInvariant, "[Error]: Use complex scalar if translational invariance is enabled");
        static_assert(!IsTransInvariant || (Dim == 1), "[Error]: Trans invariantce is not implemented in high dimension");
    private:
        ReprType repr;
        RealType dotPrec;
        FFTType planProvider;
    public:
        HubbardMatrix() = default;
        HubbardMatrix(ModelBase hubbard, ReprType repr_, RealType dotPrec_ = RealType(0));
        HubbardMatrix(const This&) = default;
        HubbardMatrix(This&&) noexcept = default;
        ~HubbardMatrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class SourceVector, class TargetVector, class Executor>
        void dot(const SourceVector& source, TargetVector& target) const;

        [[nodiscard]] ScalarType calc(size_t row, size_t col) const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using ModelBase::getHoppingT;
        using ModelBase::getRepelU;
        [[nodiscard]] const ReprType& getRepr() const noexcept { return repr; }
        [[nodiscard]] size_t getNumState() const noexcept { return repr.getNumState(); }
    protected:
        inline RealType repelElem(StateType psi) const;
        RealType hoppingElem(StateType rowPsi, StateType colPsi) const;
    private:
        template<class TargetType> void sumHopping(TargetType& target, FFTType& fft, ScalarType factor, StateType psi) const;
        template<class TargetType> void dotImpl1D(TargetType& target, ScalarType factor, size_t index) const;
        template<class TargetType> void dotImplND(TargetType& target, ScalarType factor, size_t index) const;
    };
}

namespace Physica {
    using namespace Core;

    template<class T, class U>
    class Traits<HubbardMatrix<T, U>> : public Traits<HamiltonMatrix<HubbardMatrix<T, U>>> {
    public:
        using ScalarType = T;
        using ReprType = U;

        static_assert(std::is_base_of<ReprBasis<ReprType>, ReprType>::value, "[Error]: ReprType is not a representation");
    };
}

#include "HubbardMatrixImpl.h"
#include "HubbardVecProduct.h"
