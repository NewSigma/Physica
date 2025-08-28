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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixFunction/MatrixExp.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Physics/ManyBody/Model/SquareLattice.h"

namespace Physica {
    template<Scalar T>
    class HubbardParams {
        using This = HubbardParams<T>;

        using Tr = T::RealType;
    public:
        using MatrixND = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element>;
    private:
        DenseSymmMatrix<T> hoppingMatrix;
        MatrixND expB;
        T alpha;
        T beta;
        T repelU;
        T chemMu;
        int numSplit;

        T lnCoshShift;
    public:
        HubbardParams() = default;
        template<int Dim>
        HubbardParams(T hoppingT, Tr repelU, const SquareLattice<Dim>& lattice, T beta_, T chemMu_, int numSplit_);
        HubbardParams(const This&) = default;
        HubbardParams(This&&) noexcept = default;
        ~HubbardParams() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] T calcBetaMu() const noexcept;

        [[nodiscard]] auto toDevice() const;
        [[nodiscard]] auto toDeviceAsync() const;
        void toDevice(device_obj<This>& obj) const;
        void toDeviceAsync(device_obj<This>& obj) const;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getHoppingMatrix() const noexcept { return hoppingMatrix; }
        [[nodiscard]] int getNumSite() const noexcept { return expB.getRow(); }
        [[nodiscard]] int getNumSplit() const noexcept { return numSplit; }
        [[nodiscard]] const auto& getExpB() const noexcept { return expB; }
        [[nodiscard]] T getAlpha() const noexcept { return alpha; }
        [[nodiscard]] T getBeta() const noexcept { return beta; }
        [[nodiscard]] T getRepelU() const noexcept { return repelU; }
        [[nodiscard]] T getChemMu() const noexcept { return chemMu; }
        [[nodiscard]] T getLnCoshShift() const noexcept { return lnCoshShift; }
        /* Setters */
        void setChemMu(T chemMu_);
        /* Static members */
        [[nodiscard]] static T calcAlpha(T beta, T repelU, int numSplit) noexcept;
        [[nodiscard]] static T calcBetaMu(T beta, T repelU, T chemMu) noexcept;
    private:
        template<int Dim>
        void makeHoppingMatrix(T hoppingT, const SquareLattice<Dim>& lattice);
        /* Friends */
        friend class device_obj<This>;
    };

    template<Scalar T>
    template<int Dim>
    HubbardParams<T>::HubbardParams(T hoppingT, Tr repelU_, const SquareLattice<Dim>& lattice, T beta_, T chemMu_, int numSplit_)
            : beta(beta_)
            , repelU(repelU_)
            , chemMu(chemMu_) 
            , numSplit(numSplit_) {
        static_assert(1 <= Dim && Dim <= 3, "[Error]: Invalid Dim");
        assert(!repelU.isNegative() && "[Error]: It is assumed U >= 0");
        assert(!beta.isNegative() && "[Error]: Negative temperature is invalid");
        assert(numSplit > 0 && "[Error]: Invalid NumSplit");
        makeHoppingMatrix(hoppingT, lattice);
        DenseSymmMatrix<T> hoppingMatrixB = -beta / T(numSplit) * hoppingMatrix;
        expB = exp(hoppingMatrixB);
        alpha = calcAlpha(beta, repelU, numSplit);
    }

    template<Scalar T>
    T HubbardParams<T>::calcBetaMu() const noexcept {
        return calcBetaMu(beta, repelU, chemMu);
    }

    template<Scalar T>
    void HubbardParams<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        hoppingMatrix.swap(obj.hoppingMatrix);
        expB.swap(obj.expB);
        alpha.swap(obj.alpha);
        beta.swap(obj.beta);
        repelU.swap(obj.repelU);
        chemMu.swap(obj.chemMu);
        std::swap(numSplit, obj.numSplit);

        lnCoshShift.swap(obj.lnCoshShift);
    }

    template<Scalar T>
    void HubbardParams<T>::setChemMu(T chemMu_) {
        chemMu = std::move(chemMu_);
        lnCoshShift = lncosh(calcBetaMu());
    }

    template<Scalar T>
    T HubbardParams<T>::calcAlpha(T beta, T repelU, int numSplit) noexcept {
        const T betaM = beta / T(numSplit);
        const T x = betaM * repelU;
        return x * T(0.5) + ln1p(sqrt(T(1) - exp(-x)));
    }

    template<Scalar T>
    T HubbardParams<T>::calcBetaMu(T beta, T repelU, T chemMu) noexcept {
        return beta * (chemMu - repelU * T(0.5));
    }

    template<Scalar T>
    template<int Dim>
    void HubbardParams<T>::makeHoppingMatrix(T hoppingT, const SquareLattice<Dim>& lattice) {
        const size_t numSite = lattice.getNumSuperCellSite();
        hoppingMatrix.resize(numSite);
        hoppingMatrix = T(0);

        for (size_t from = 0; from < numSite; ++from) {
            if constexpr (SquareLattice<Dim>::UntrivialNearestNeighbor) {
                const auto& targets = lattice.getHopIndexArray()[from];
                for (size_t to : targets)
                    hoppingMatrix(from, to) -= hoppingT;
            }
            else
                hoppingMatrix(from, (from + 1) % numSite) = -hoppingT;
        }
    }
}
