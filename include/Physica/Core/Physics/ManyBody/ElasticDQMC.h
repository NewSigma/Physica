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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiagMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.h"
#include "DQMCImpl/ImagKinetic.h"
#include "FreqDQMC.h"

namespace Physica {
    template<Scalar T>
    class ElasticDQMC {
        using This = ElasticDQMC<T>;
        using Params = HubbardParams<T>;
        using GreenPair = ImagKinetic<T>::GreenPair;

        using Tr = T::RealType;
        using Tc = T::ComplexType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;

        constexpr static bool isComplex = T::isComplex;
    private:
        const Params* params;
        VectorND<Trv> rSquareOmegas;
        GreenPair greens;
        GreenPair buffer;

        MatrixND<T> actionR;
        SymmEigenSolver<T> eig;

        T lnAbsDet;
    public:
        ElasticDQMC() = delete;
        ElasticDQMC(const Params& params_, Trv freqDensity);
        ElasticDQMC(const This&) = default;
        ElasticDQMC(This&&) noexcept = default;
        ~ElasticDQMC() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<RNG R>
        void step_random();

        template<RNG R>
        void step();
        template<RNG R>
        void step_for(int numStep);

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getParams() const noexcept { return *params; }
        [[nodiscard]] int getNumSite() const noexcept { return params->getNumSite(); }
        [[nodiscard]] int getNumSplit() const noexcept { return params->getNumSplit(); }
        [[nodiscard]] auto& getGreens() noexcept { return greens; }
        [[nodiscard]] constexpr static Trv getSign() noexcept { return 1; }
        [[nodiscard]] constexpr static Trv getRSign() noexcept { return 1; }
    private:
        template<RNG R>
        [[nodiscard]] T randAuxField();
        [[nodiscard]] T diagonalize();
    };

    template<Scalar T>
    ElasticDQMC<T>::ElasticDQMC(const Params& params_, Trv freqDensity)
            : params(&params_)
            , rSquareOmegas(FreqDQMC<Tc>::calcFreqCutoff(params_.getBeta(), freqDensity))
            , greens(2, params_.getNumSite())
            , buffer(2, params_.getNumSite())
            , actionR(params_.getHoppingMatrix() * params_.getBeta())
            , eig(params_.getNumSite(), true) {
        assert(freqDensity.isPositive());
        for (size_t i = 0; i < rSquareOmegas.getLength(); ++i)
            rSquareOmegas[i] = T(2 * i + 1);
        rSquareOmegas *= MathConst<T>::pi;
        rSquareOmegas = reciprocal(square(rSquareOmegas));
    }

    template<Scalar T>
    template<RNG R>
    void ElasticDQMC<T>::step_random() {
        for (int i = 0; i < getNumSite(); ++i)
            actionR(i, i) = randAuxField<R>();
        lnAbsDet = diagonalize();
        buffer.swap(greens);
    }

    template<Scalar T>
    template<RNG R>
    void ElasticDQMC<T>::step() {
        const int site = std::uniform_int_distribution<int>(0, getNumSite() - 1)(R::getInstance());
        const T save = std::exchange(actionR(site, site), randAuxField<R>());

        const auto lnAbsDet1 = diagonalize();
        const bool accept = Trv::template random_uniform<R>() < exp(lnAbsDet1 - lnAbsDet);
        if (accept) {
            lnAbsDet = lnAbsDet1;
            buffer.swap(greens);
        }
        else
            actionR(site, site) = save;
    }

    template<Scalar T>
    template<RNG R>
    void ElasticDQMC<T>::step_for(int numStep) {
        assert(numStep >= 0 && "[Error]: Invalid step num");
        for (int _ = 0; _ < numStep; ++_)
            step<R>();
    }

    template<Scalar T>
    void ElasticDQMC<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(params, obj.params);
        rSquareOmegas.swap(obj.rSquareOmegas);
        greens.swap(obj.greens);

        actionR.swap(obj.actionR);
        eig.swap(obj.eig);

        lnAbsDet.swap(obj.lnAbsDet);
    }

    template<Scalar T>
    template<RNG R>
    T ElasticDQMC<T>::randAuxField() {
        Tr betaU = params->getBeta() * params->getRepelU();
        Vector3D<Tr> shifts{-betaU, 0, betaU};
        const Tr shift = shifts[std::uniform_int_distribution<int>(0, 2)(R::getInstance())];
        return T::template random_normal<R>() * sqrt(betaU) + shift;
    }

    template<Scalar T>
    auto ElasticDQMC<T>::diagonalize() -> T {
        T lnAbsDet = 0;
        for (int spin : {0, 1}) {
            const Tr shift = params->getBeta() * fma(params->getRepelU(), Tr(-0.5), params->getChemMu());
            using Tf = Diff<T, DiffMode::Reverse>;
            eig.compute(actionR);
            auto eigenvalues = VectorND<Tf>(eig.getEigenvalues());
            lnAbsDet += ln1p_elem(square(eigenvalues - shift) * rSquareOmegas.transpose()).sum().reverse();
            actionR.diag() = -actionR.diag();

            auto& green = buffer[spin];
            MatrixND<T> temp = eig.getEigenvectors() * DiagMatrix<T>(eigenvalues.grads());
            green = temp * eig.getEigenvectors().transpose();
            green.diag() = Trv(0.5) + green.diag();
        }
        Tr betaU = params->getBeta() * params->getRepelU();
        lnAbsDet -= ln1pexp(lncosh(actionR.diag()) + fma(betaU, Trv(-0.5), MathConst<Trv>::ln2)).sum();
        return lnAbsDet;
    }
}
