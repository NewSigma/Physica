/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Transform/FFT.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/ReprBasis.h"
#include "Physica/Core/Physics/ManyBody/Model/SquareLattice.h"
#include "HamiltonMatrix.h"

namespace Physica {
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
        using typename Lattice::IndexType;
        using typename Lattice::ArgVector;
    public:
        using typename Base::StateType;
    protected:
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
        FFT1D plan;
    public:
        HubbardMatrix() = default;
        HubbardMatrix(T hoppingT_, Tr repelU_, Lattice lattice, Repr repr_);
        HubbardMatrix(const This&) = default;
        HubbardMatrix(This&&) noexcept = default;
        ~HubbardMatrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }

        [[nodiscard]] auto operator*(this auto&&, Vector auto&& v) noexcept;
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;
        [[nodiscard]] T trace() const;
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Lattice::getNeighbors;
        [[nodiscard]] T getHoppingT() const noexcept { return hoppingT; }
        [[nodiscard]] Tr getRepelU() const noexcept { return repelU; }
        [[nodiscard]] const Lattice& getLattice() const noexcept { return *this; }
        [[nodiscard]] const auto& getRepr() const noexcept { return repr; }
        [[nodiscard]] const auto& getPlan() const noexcept { return plan; }
        /* Setters */
        void setPhaseArgs(ArgVector phaseArgs) noexcept requires(BC == BoundaryCond::TBC);
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return MatrixMajor::BothMajor; }
    protected:
        Tr repelElem(StateType psi) const noexcept;
        T hoppingElem(StateType rowPsi, StateType colPsi) const noexcept;
    private:
        Vector2D<T> calcBoundaryPhase(int site, int site1) const noexcept;
        template<Matrix, Vector> friend class GEMV;
    };
}

namespace Physica {
    template<Scalar T, Representation Repr, BoundaryCond BC>
    class Traits<HubbardMatrix<T, Repr, BC>> {
    public:
        using ScalarType = T;
        using ReprType = Repr;
        constexpr static BoundaryCond Boundary = BC;
    };
}

#include "HubbardMatrixImpl.h"
#include "HubbardGEMV.h"
