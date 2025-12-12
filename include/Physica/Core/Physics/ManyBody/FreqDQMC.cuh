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

#include "FreqDQMC.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/DenseLU.cuh"

namespace Physica {
    template<Scalar T>
    class device_obj<FreqDQMC<T>> {
        using host_obj = FreqDQMC<T>;
        using This = device_obj<host_obj>;
        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;

        using GreenPair = ImagKinetic<Tr>::GreenPair;
    private:
        Array<device_obj<DenseLU<T, false>>, 2> lu;
        ActionMatrix<T> action;
        device_obj<MatrixND<T>> linearRHS;
        MatrixND<T> detBuffer;
        MatrixND<T> solBuffer;

        GreenPair greens;
        Trv lnAbsDet = Trv::nan();
        Trv sign = 1;
    public:
        device_obj(const HubbardParams<Tr>& params, Trv freqDensity);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<RNG R>
        void step_random();

        template<RNG R>
        void step();
        template<RNG R>
        void step_for(int numStep);
        /* Getters */
        [[nodiscard]] const auto& getParams() const noexcept { return action.getParams(); }
        [[nodiscard]] int getNumSite() const noexcept { return getParams().getNumSite(); }
        [[nodiscard]] int getNumFreq() const noexcept { return action.getNumFreq(); }
        [[nodiscard]] int getNumSplit() const noexcept { return getParams().getNumSplit(); }
        [[nodiscard]] auto& getGreens() noexcept { return greens; }
        [[nodiscard]] Trv getSign() const noexcept { return sign; }
        [[nodiscard]] Trv getRSign() const noexcept { return getSign(); }
        /* Static members */
        [[nodiscard]] static int calcFreqCutoff(Trv beta, Trv freqDensity) noexcept;
    private:
        [[nodiscard]] Vector2D<Trv> calcDet();
        void calcGreen();
    };

    template<Scalar T>
    device_obj<FreqDQMC<T>>::device_obj(const HubbardParams<Tr>& params, Trv freqDensity)
            : action(params, calcFreqCutoff(params.getBeta(), freqDensity))
            , greens(2, params.getNumSite()) {
        size_t size = action.getOrder();
        for (auto& spinLU : lu)
            spinLU.resize(size);
        linearRHS.resize(size);
        detBuffer.resize(size);
        solBuffer.resize(size);

        action.assign_kinetic(detBuffer);
    }

    template<Scalar T>
    template<RNG R>
    void device_obj<FreqDQMC<T>>::step_random() {
        action.template randAuxField<R>();
        auto [lnAD, sgnD] = calcDet();
        lnAbsDet = lnAD;
        sign = sgnD;
        calcGreen();
    }

    template<Scalar T>
    template<RNG R>
    void device_obj<FreqDQMC<T>>::step() {
        assert(lnAbsDet.isFinite() && "[Error]: Should random initialize before monte carlo step");
        const int site = std::uniform_int_distribution<int>(0, getNumSite() - 1)(R::getInstance());
        auto aux = action.getAuxField().col(site);
        const VectorND<T> save = aux;
        action.template randAuxField<R>(site);

        auto [lnAD, sgnD] = calcDet();
        const bool accept = Trv::template random_uniform<R>() < exp(lnAD - lnAbsDet);
        if (accept) {
            lnAbsDet = lnAD;
            sign = sgnD;
            calcGreen();
        }
        else
            aux = save;
    }

    template<Scalar T>
    template<RNG R>
    void device_obj<FreqDQMC<T>>::step_for(int numStep) {
        for (int i = 0; i < numStep; ++i)
            step<R>();
    }

    template<Scalar T>
    int device_obj<FreqDQMC<T>>::calcFreqCutoff(Trv beta, Trv freqDensity) noexcept {
        assert(!freqDensity.isNegative());
        int i = static_cast<int>((beta * freqDensity).toMachine() * 0.25);
        return (i + 1) * 4; // Multiple of 4 so that SIMD works
    }

    template<Scalar T>
    auto device_obj<FreqDQMC<T>>::calcDet() -> Vector2D<Trv> {
        action.assign_potential(detBuffer);
        lu[0].compute(detBuffer);

        action.getAuxField() = -action.getAuxField();
        action.assign_potential(detBuffer);
        lu[1].compute(detBuffer);

        Trv lnAD = 0, sgnD = 1;
        for (const auto& spinLU : lu) {
            lnAD += spinLU.lnAbsDet();
            sgnD *= spinLU.sgndet().real();
        }
        return {lnAD, sgnD};
    }

    template<Scalar T>
    void device_obj<FreqDQMC<T>>::calcGreen() {
        const int numSite = getNumSite();
        for (int spin : {0, 1}) {
            auto& spinLU = lu[spin];
            solBuffer = UnitMatrix<T>(action.getOrder());
            solBuffer.toDevice(linearRHS);
            spinLU.solve(linearRHS);
            linearRHS.toHost(solBuffer);

            auto& green = greens[spin];
            green.zeros();
            for (int _ = 0, offset = 0; _ < action.getNumFreq() * 2; ++_) {
                green += solBuffer.block(offset, numSite, offset, numSite).reals();
                offset += numSite;
            }
            green *= reciprocal(getParams().getBeta());
            green.diag() += Trv(0.5);
        }
    }
}
