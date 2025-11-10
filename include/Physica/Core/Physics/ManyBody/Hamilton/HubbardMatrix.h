/*
 * Copyright 2024-2025 Weibo He.
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
#include "Physica/Core/Physics/ManyBody/Model/SquareLattice.h"
#include "HamiltonMatrix.h"

namespace Physica {
    template<Scalar, Representation, BoundaryCond, Vector>
    class HubbardVecProd;
    /**
     * Refer to [1] for applied symmetries
     * 
     * Reference:
     * [1] J. Korean Phys. Soc. 76, 670–683 (2020); https://doi.org/10.3938/jkps.76.670
     */
    template<Scalar T, Representation Repr, BoundaryCond BC = BoundaryCond::PBC>
    class HubbardMatrix : public HamiltonMatrix<HubbardMatrix<T, Repr, BC>>
                        , public SquareLattice<Repr::Dim, BC> {
        using This = HubbardMatrix<T, Repr, BC>;
        using Base = HamiltonMatrix<This>;
        using Lattice = SquareLattice<Repr::Dim, BC>;

        using typename Base::StateType;
        using typename Lattice::IndexType;
        using typename Lattice::ArgVector;

        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Tc = T::ComplexType;
    public:
        using Base::Dim;
        using Base::NumSite;
        using FFT1D = FFT<Tr, 1>;
    private:
        T hoppingT;
        Tr repelU;
        Repr repr;
        DenseVector<Tc, Dim> phases;
        FFT1D planProvider;
    public:
        HubbardMatrix() = default;
        HubbardMatrix(T hoppingT_, Tr repelU_, Lattice lattice, Repr repr_);
        HubbardMatrix(const This&) = default;
        HubbardMatrix(This&&) noexcept = default;
        ~HubbardMatrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }

        [[nodiscard]] auto operator*(const Vector auto& v) const noexcept;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;
        [[nodiscard]] Tv calc_value(size_t row, size_t col) const;
        [[nodiscard]] T trace() const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Lattice::getNeighbors;
        [[nodiscard]] T getHoppingT() const noexcept { return hoppingT; }
        [[nodiscard]] Tr getRepelU() const noexcept { return repelU; }
        [[nodiscard]] const Lattice& getLattice() const noexcept { return *this; }
        [[nodiscard]] const auto& getRepr() const noexcept { return repr; }
        /* Setters */
        void setPhaseArgs(ArgVector phaseArgs) noexcept requires(BC == BoundaryCond::TBC);
    protected:
        Tr repelElem(StateType psi) const noexcept;
        T hoppingElem(StateType rowPsi, StateType colPsi) const noexcept;
    private:
        Vector2D<T> calcBoundaryPhase(int site, int site1) const noexcept;
        template<Scalar, Representation, BoundaryCond, Vector>
        friend class HubbardVecProd;
    };
}

namespace Physica {
    template<Scalar T, Representation Repr, BoundaryCond BC>
    class Traits<HubbardMatrix<T, Repr, BC>> : public Traits<HamiltonMatrix<HubbardMatrix<T, Repr, BC>>> {
    public:
        using ScalarType = T;
        using ReprType = Repr;
        constexpr static BoundaryCond Boundary = BC;
    };
}

#include "HubbardMatrixImpl.h"
#include "HubbardVecProd.h"
