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

#include "Physica/Core/Math/Transform/FFT.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/ReprBasis.h"
#include "Physica/Core/Physics/ManyBody/Model/Hubbard.h"
#include "HamiltonMatrix.h"

namespace Physica::Core {
    class ThreadExecutor;
    /**
     * Refer to [1] for applied symmetries
     * 
     * Reference:
     * [1] J. Korean Phys. Soc. 76, 670–683 (2020); https://doi.org/10.3938/jkps.76.670
     */
    template<Scalar T, Representation U>
    class HubbardMatrix
            : public HamiltonMatrix<HubbardMatrix<T, U>>
            , public Hubbard<typename T::RealType, U::Dim> {
        using RealType = T::RealType;
        using This = HubbardMatrix<T, U>;
        using Base = HamiltonMatrix<This>;
        using ModelBase = Hubbard<RealType, U::Dim>;
        
        using StateType = U::StateType;
        using typename ModelBase::IndexType;
    public:
        using FFTType = FFT<RealType, 1>;
        constexpr static unsigned int Dim = U::Dim;
        constexpr static unsigned int NumSite = StateType::NumSite;
        constexpr static unsigned int SiteDOF = StateType::SiteDOF;
        constexpr static bool IsTransInvariant = Traits<U>::IsTransInvariant;
        static_assert((IsTransInvariant && T::isComplex) || !IsTransInvariant, "[Error]: Use complex scalar if translational invariance is enabled");
        static_assert(!IsTransInvariant || (Dim == 1), "[Error]: Trans invariantce is not implemented in high dimension");
    private:
        U repr;
        FFTType planProvider;
    public:
        HubbardMatrix() = default;
        HubbardMatrix(ModelBase hubbard, U repr_);
        HubbardMatrix(const This&) = default;
        HubbardMatrix(This&&) noexcept = default;
        ~HubbardMatrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class Functor>
        void forNeighSites(Functor func, int site) const;

        [[nodiscard]] T calc(size_t row, size_t col) const;
        [[nodiscard]] T trace() const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using ModelBase::getHoppingT;
        using ModelBase::getRepelU;
        using ModelBase::getHopIndexArray;
        [[nodiscard]] const ModelBase& getModel() const noexcept { return *this; }
        [[nodiscard]] const U& getRepr() const noexcept { return repr; }
    protected:
        inline RealType repelElem(StateType psi) const;
        RealType hoppingElem(StateType rowPsi, StateType colPsi) const;
    private:
        template<Matrix, Vector> friend class MatrixVectorProduct;
    };
}

namespace Physica {
    template<Scalar T, Representation U>
    class Traits<HubbardMatrix<T, U>> : public Traits<HamiltonMatrix<HubbardMatrix<T, U>>> {
    public:
        using ScalarType = T;
        using ReprType = U;
    };
}

#include "HubbardMatrixImpl.h"
#include "HubbardVecProduct.h"
