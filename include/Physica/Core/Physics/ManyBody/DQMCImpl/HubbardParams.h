/*
 * Copyright 2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixFunction/MatrixExp.h"
#include "Physica/Core/Physics/ManyBody/Hamilton/HubbardMatrix.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/FermiRepr.h"

namespace Physica {
    template<Scalar T>
    class HubbardParams {
        using This = HubbardParams<T>;

        using Tr = T::RealType;
        using Tc = T::ComplexType;
    public:
        using MatrixND = DenseMatrix<T>;
    private:
        MatrixND hoppingMatrix;
        MatrixND expT;
        Tr alpha;
        Tr beta;
        Tr repelU;
        Tr chemMu;
        int numSplit = 0;

        VectorND<Tr> lnSpinWeights;
    public:
        HubbardParams() = default;
        template<int Dim, BoundaryCond BC>
        HubbardParams(Tr hoppingT, Tr repelU_, const SquareLattice<Dim, BC>& lattice, Tr beta_, Tr chemMu_, int numSplit_);
        HubbardParams(const This&) = default;
        HubbardParams(This&&) noexcept = default;
        ~HubbardParams() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] Tr calcBetaMu() const noexcept;
        [[nodiscard]] MatrixND calcInvExpT() const;

        [[nodiscard]] auto toDevice() const;
        [[nodiscard]] auto toDeviceAsync() const;
        void toDevice(device_obj<This>& obj) const;
        void toDeviceAsync(device_obj<This>& obj) const;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getHoppingMatrix() const noexcept { return hoppingMatrix; }
        [[nodiscard]] int getNumSite() const noexcept { return expT.getRow(); }
        [[nodiscard]] int getNumSplit() const noexcept { return numSplit; }
        [[nodiscard]] const auto& getExpT() const noexcept { return expT; }
        [[nodiscard]] Tr getAlpha() const noexcept { return alpha; }
        [[nodiscard]] Tr getBeta() const noexcept { return beta; }
        [[nodiscard]] Tr getRepelU() const noexcept { return repelU; }
        [[nodiscard]] Tr getChemMu() const noexcept { return chemMu; }
        [[nodiscard]] const auto& getLnSpinWeights() const noexcept { return lnSpinWeights; }
        /* Static members */
        [[nodiscard]] static Tr calcAlpha(Tr beta, Tr repelU, int numSplit) noexcept;
        [[nodiscard]] static Tr calcBetaMu(Tr beta, Tr repelU, Tr chemMu) noexcept;
    private:
        template<int Dim, BoundaryCond BC>
        void makeHoppingMatrix(Tr hoppingT, const SquareLattice<Dim, BC>& lattice);
        void makeLnSpinWeights();
        /* Friends */
        friend class device_obj<This>;
    };

    template<Scalar T>
    template<int Dim, BoundaryCond BC>
    HubbardParams<T>::HubbardParams(Tr hoppingT, Tr repelU_, const SquareLattice<Dim, BC>& lattice, Tr beta_, Tr chemMu_, int numSplit_)
            : beta(beta_)
            , repelU(repelU_)
            , chemMu(chemMu_)
            , numSplit(numSplit_)
            , lnSpinWeights(lattice.getNumSuperCellSite()) {
        assert(!repelU.isNegative() && "[Error]: It is assumed U >= 0");
        assert(!beta.isNegative() && "[Error]: Negative temperature is invalid");
        assert(numSplit > 0 && "[Error]: Invalid NumSplit");
        makeHoppingMatrix<Dim, BC>(hoppingT, lattice);

        MatrixND hoppingMatrixB = -beta / T(numSplit) * hoppingMatrix;
        expT = exp(hoppingMatrixB);
        alpha = calcAlpha(beta, repelU, numSplit);
        makeLnSpinWeights();
    }

    template<Scalar T>
    auto HubbardParams<T>::calcBetaMu() const noexcept -> Tr {
        return calcBetaMu(beta, repelU, chemMu);
    }

    template<Scalar T>
    auto HubbardParams<T>::calcInvExpT() const -> MatrixND {
        MatrixND hoppingMatrixB = beta / T(numSplit) * hoppingMatrix;
        return exp(hoppingMatrixB);
    }

    template<Scalar T>
    void HubbardParams<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        hoppingMatrix.swap(obj.hoppingMatrix);
        expT.swap(obj.expT);
        alpha.swap(obj.alpha);
        beta.swap(obj.beta);
        repelU.swap(obj.repelU);
        chemMu.swap(obj.chemMu);
        std::swap(numSplit, obj.numSplit);
    }

    template<Scalar T>
    auto HubbardParams<T>::calcAlpha(Tr beta, Tr repelU, int numSplit) noexcept -> Tr {
        const Tr betaM = beta / Tr(numSplit);
        const Tr x = betaM * repelU;
        return x * Tr(0.5) + ln1p(sqrt(Tr(1) - exp(-x)));
    }

    template<Scalar T>
    auto HubbardParams<T>::calcBetaMu(Tr beta, Tr repelU, Tr chemMu) noexcept -> Tr {
        return beta * (chemMu - repelU * Tr(0.5));
    }

    template<Scalar T>
    template<int Dim, BoundaryCond BC>
    void HubbardParams<T>::makeHoppingMatrix(Tr hoppingT, const SquareLattice<Dim, BC>& lattice) {
        const auto CheckTemplateParam = HubbardMatrix<T, FermiRepr<Dim, Dynamic, false>, BC>{};

        const size_t numSite = lattice.getNumSuperCellSite();
        hoppingMatrix.resize(numSite);
        hoppingMatrix.zeros();

        using enum BoundaryCond;
        if constexpr (BC == PBC) {
            for (size_t from = 0; from < numSite; ++from) {
                if constexpr (Dim > 1) {
                    const auto& targets = lattice.getNeighbors()[from];
                    for (size_t to : targets) {
                        hoppingMatrix(from, to) -= hoppingT;
                        hoppingMatrix(to, from) -= hoppingT;
                    }
                }
                else {
                    size_t next = (from + 1) % numSite;
                    hoppingMatrix(from, next) = -hoppingT;
                    hoppingMatrix(next, from) = -hoppingT;
                }
            }
        }
        else {
            static_assert(BC == APBC || BC == TBC, "[Error]: Not implemented");
            const auto phases = lattice.template calcPhase<Tr>();
            using PhaseScalar = decltype(phases)::ScalarType;

            const auto& boundary = lattice.getSiteBoundaryMap();
            for (size_t from = 0; from < numSite; ++from) {
                if constexpr (Dim > 1) {
                    const auto& targets = lattice.getNeighbors()[from];
                    for (size_t to : targets) {
                        PhaseScalar phase = 1;
                        auto pair = std::make_pair(from, to);
                        if (boundary.contains(pair)) {
                            int dim = boundary.find(pair)->second;
                            phase = phases[dim];
                        }
                        hoppingMatrix(from, to) -= hoppingT * phase;
                        hoppingMatrix(to, from) -= hoppingT * phase.conjugate();
                    }
                }
                else {
                    size_t next = (from + 1) % numSite;
                    PhaseScalar phase = 1;
                    auto pair = std::make_pair(from, next);
                    if (boundary.contains(pair)) {
                        int dim = boundary.find(pair)->second;
                        phase = phases[dim];
                    }
                    hoppingMatrix(from, next) = -hoppingT * phase;
                    hoppingMatrix(next, from) = -hoppingT * phase.conjugate();
                }
            }
        }
    }

    template<Scalar T>
    void HubbardParams<T>::makeLnSpinWeights() {
        const Tr lncoshBetaMu = lncosh(calcBetaMu());
        const int numSplit = getNumSplit();
        lnSpinWeights.resize(numSplit + 1);
        for (int spinUp = 0; spinUp <= numSplit; ++spinUp) {
            int sumSpin = 2 * spinUp - numSplit;
            lnSpinWeights[spinUp] = ln1pexp(lncosh(alpha * sumSpin) - lncoshBetaMu);
        }
    }
}
