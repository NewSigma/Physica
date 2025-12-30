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

#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/OpenModel.h"

namespace Physica {
    template<Scalar T>
    class HamiltonMC {
        using This = HamiltonMC<T>;
        using Engine = RPMD<T, 1, 1>;
        using KineticModel = OpenModel<T, 1, 1>;
        using Tr = T::RealType;
        using Trv = Tr::ValueType;

        Engine engine;
        KineticModel kinetic;
        T duration;
        VectorND<T> prevX;

        int64_t numAccept = 0;
        int64_t numTotal = 0;
    public:
        HamiltonMC(VectorND<T> init, T stepsize, T duration, T temperatureT = 1);
        HamiltonMC(const This&) = default;
        HamiltonMC(This&&) noexcept = default;
        ~HamiltonMC() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        template<RNG R>
        [[nodiscard]] const auto& step(auto&& forceModel);
        /* Getters */
        [[nodiscard]] const auto& getEngine() const noexcept { return engine; }
        [[nodiscard]] size_t getDOF() const noexcept { return engine.getDOF(); }
        [[nodiscard]] Trv getAcceptRate() const noexcept { return Trv(numAccept) / Trv(numTotal); }
    private:
        static MDCell<T, 1> makeDummyCell(size_t dof);
    };

    template<Scalar T>
    HamiltonMC<T>::HamiltonMC(VectorND<T> init, T stepsize, T duration, T temperatureT)
            : engine(makeDummyCell(init.getLength()), 1, 1, temperatureT, stepsize)
            , kinetic(temperatureT, 1)
            , duration(duration)
            , prevX(std::move(init)) {
        auto col = engine.getPhaseMatrix().col(0);
        col.head(getDOF()).zeros();
        col.tail(getDOF()) = prevX;
    }

    template<Scalar T>
    template<RNG R>
    const auto& HamiltonMC<T>::step(auto&& forceModel) {
        engine.template initMomentum<KineticModel, R>();

        const T prevE = engine.template calcKinetic<KineticModel>() + engine.calcPotential(forceModel);
        engine.nve_step_for(duration, kinetic, forceModel);

        const T nowE = engine.template calcKinetic<KineticModel>() + engine.calcPotential(forceModel);
        const bool accept = T::template random_uniform<R>() < exp((prevE - nowE) / engine.getTemperature());
        if (accept) [[likely]] {
            prevX = engine.getPhaseMatrix().col(0).tail(getDOF());
            numAccept += 1;
        }
        else
            engine.getPhaseMatrix().col(0).tail(getDOF()) = prevX;

        numTotal += 1;
        return prevX;
    }

    template<Scalar T>
    MDCell<T, 1> HamiltonMC<T>::makeDummyCell(size_t dof) {
        return MDCell<T, 1>(dof, VectorND<T>(dof, 1));
    }
}
