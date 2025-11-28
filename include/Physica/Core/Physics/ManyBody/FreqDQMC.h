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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/SymmEigenSolver.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.h"
#include "Physica/Core/Scalar/Diff.h"
#include "DQMCImpl/ImagKinetic.h"
#include "DQMCImpl/GreenProd.h"

namespace Physica {
    template<Scalar T>
    class FreqDQMC {
        using This = FreqDQMC<T>;
        using Params = HubbardParams<T>;

        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;

        constexpr static bool isComplex = T::isComplex;
    private:
        const Params* params;
        VectorND<Trv> rSquareOmegas;
        ImagKinetic<T> kinetic;

        MatrixND<T> actionR;
        SymmEigenSolver<T> eig;

        T lnZ;
    public:
        FreqDQMC() = delete;
        FreqDQMC(const Params& params_, Trv freqDensity);
        FreqDQMC(const This&) = default;
        FreqDQMC(This&&) noexcept = default;
        ~FreqDQMC() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<RNG R = Random<>>
        void step_random();

        template<RNG R = Random<>>
        void step_spin();
        template<RNG R = Random<>>
        void step_spin_for(int numStep);

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getParams() const noexcept { return *params; }
        [[nodiscard]] int getNumSite() const noexcept { return params->getNumSite(); }
        [[nodiscard]] int getNumSplit() const noexcept { return params->getNumSplit(); }
        [[nodiscard]] auto& getGreens() noexcept { return kinetic.getGreens(); }
        /* Static members */
        [[nodiscard]] static Trv getSign() noexcept { return 1; }
    private:
        /* Operations */
        void metropolis(int site, Tr prob);
        template<RNG R>
        [[nodiscard]] T randAuxField();
        [[nodiscard]] T calcLnZ();
        /* Static members */
        [[nodiscard]] static int calcFreqCutoff(Trv beta, Trv freqDensity);
    };

    template<Scalar T>
    FreqDQMC<T>::FreqDQMC(const Params& params_, Trv freqDensity)
            : params(&params_)
            , rSquareOmegas(calcFreqCutoff(params_.getBeta(), freqDensity))
            , kinetic(params_.getNumSite(), params_.getNumSplit())
            , actionR(params_.getHoppingMatrix())
            , eig(params_.getNumSite(), true) {
        assert(freqDensity.isPositive());
        for (size_t i = 0; i < rSquareOmegas.getLength(); ++i)
            rSquareOmegas[i] = T(2 * i + 1);
        rSquareOmegas *= MathConst<T>::pi / params->getBeta();
        rSquareOmegas = reciprocal(square(rSquareOmegas));
    }

    template<Scalar T>
    template<RNG R>
    void FreqDQMC<T>::step_random() {
        for (int i = 0; i < getNumSite(); ++i)
            actionR(i, i) = randAuxField<R>();
        lnZ = calcLnZ().value();
    }

    template<Scalar T>
    template<RNG R>
    void FreqDQMC<T>::step_spin() {
        const int site = std::uniform_int_distribution<int>(0, getNumSite() - 1)(R::getInstance());
        const T save = std::exchange(actionR(site, site), randAuxField<R>());

        const auto lnZ1 = calcLnZ();
        const bool accept = Trv::template random_uniform<R>() < exp(lnZ1.value() - lnZ.value());
        if (accept)
            lnZ = lnZ1;
        else
            actionR(site, site) = save;
    }

    template<Scalar T>
    template<RNG R>
    void FreqDQMC<T>::step_spin_for(int numStep) {
        assert(numStep >= 0 && "[Error]: Invalid step num");
        for (int _ = 0; _ < numStep; ++_)
            step_spin<R>();
    }

    template<Scalar T>
    void FreqDQMC<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(params, obj.params);
        rSquareOmegas.swap(obj.rSquareOmegas);
        kinetic.swap(obj.kinetic);

        actionR.swap(obj.actionR);
        eig.swap(obj.eig);

        lnZ.swap(obj.lnZ);
    }

    template<Scalar T>
    template<RNG R>
    T FreqDQMC<T>::randAuxField() {
        return T::template random_normal<R>() * sqrt(params->getRepelU() / params->getBeta());
    }

    template<Scalar T>
    auto FreqDQMC<T>::calcLnZ() -> T {
        auto diagonalize = [&](MatrixND<T>& green) {
            using Tf = Diff<T, DiffMode::Reverse>;
            const T shiftMu = params->getChemMu() - params->getRepelU() * 0.5;
            const T repBeta = reciprocal(params->getBeta());

            eig.compute(actionR);
            auto eigenvalues = VectorND<Tf>(eig.getEigenvalues());
            T lnZ = ln1p_elem(square(eigenvalues - shiftMu) * rSquareOmegas.transpose()).sum().reverse(repBeta);
            MatrixND<T> temp = eig.getEigenvectors() * DiagMatrix<T>(eigenvalues.grads());
            green = temp * eig.getEigenvectors().transpose();
            green.diag() = Trv(0.5) + green.diag();
            return lnZ;
        };

        T lnZU = diagonalize(getGreens()[0]);

        actionR.diag() = -actionR.diag();
        T lnZD = diagonalize(getGreens()[1]);
        return lnZU + lnZD;
    }

    template<Scalar T>
    int FreqDQMC<T>::calcFreqCutoff(Trv beta, Trv freqDensity) {
        int i = static_cast<int>((beta * freqDensity).toMachine() * 0.25);
        return (i + 1) * 4; // Multiple of 4 so that SIMD works
    }
}
