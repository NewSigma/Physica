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

namespace Physica {
    template<Scalar, Representation, Vector>
    class HubbardVecProd;
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
        using This = HubbardMatrix<T, U>;
        using Base = HamiltonMatrix<This>;
        using Tr = T::RealType;
        using ModelBase = Hubbard<Tr, U::Dim>;
        
        using typename Base::StateType;
        using typename ModelBase::IndexType;
    public:
        using FFTType = FFT<Tr, 1>;
        using Base::NumSite;
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

        template<Vector V>
        [[nodiscard]] auto operator*(const V& v) const noexcept;
        /* Operations */
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
        Tr repelElem(StateType psi) const;
        Tr hoppingElem(StateType rowPsi, StateType colPsi) const;
    private:
        template<Scalar, Representation, Vector>
        friend class HubbardVecProd;
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
#include "HubbardVecProd.h"
