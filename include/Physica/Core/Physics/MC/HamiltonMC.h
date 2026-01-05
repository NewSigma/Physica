/*
 * Copyright 2025-2026 Weibo He.
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
        HamiltonMC(VectorND<T> mass, T stepsize, T duration);
        HamiltonMC(const This&) = default;
        HamiltonMC(This&&) noexcept = default;
        ~HamiltonMC() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        /* Operations */
        template<RNG R, ExecutePolicy P = Sequential>
        [[nodiscard]] const auto& step(auto&& forceModel);

        void reset() noexcept;
        /* Getters */
        [[nodiscard]] auto& getEngine() noexcept { return engine; }
        [[nodiscard]] size_t getDOF() const noexcept { return engine.getDOF(); }
        [[nodiscard]] const auto& getPrevX() const noexcept { return prevX; }
        [[nodiscard]] Trv getAcceptRate() const noexcept { return Trv(numAccept) / Trv(numTotal); }
        /* Setters */
        void setInitial(VectorND<T> init);
    private:
        static MDCell<T, 1> makeDummyCell(VectorND<T> mass);
    };

    template<Scalar T>
    HamiltonMC<T>::HamiltonMC(VectorND<T> mass, T stepsize, T duration)
            : engine(makeDummyCell(std::move(mass)), 1, 1, 1, stepsize)
            , kinetic(1, 1)
            , duration(duration) {
        assert(duration > stepsize);
        engine.getPhaseMatrix().zeros();
    }

    template<Scalar T>
    template<RNG R, ExecutePolicy P>
    const auto& HamiltonMC<T>::step(auto&& forceModel) {
        engine.template initMomentum<KineticModel, R>();

        const T prevE = engine.template calcClassicalInternalEnergy<P>(forceModel);
        engine.template nve_step_for<P>(duration, kinetic, forceModel);

        const T curE = engine.template calcClassicalInternalEnergy<P>(forceModel);
        const bool accept = T::template random_uniform<R>() < exp((prevE - curE) / engine.getTemperature());
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
    void HamiltonMC<T>::reset() noexcept {
        numAccept = 0;
        numTotal = 0;
    }

    template<Scalar T>
    void HamiltonMC<T>::setInitial(VectorND<T> init) {
        assert(getDOF() == init.getLength());
        prevX = std::move(init);

        auto col = engine.getPhaseMatrix().col(0);
        col.head(getDOF()).zeros();
        col.tail(getDOF()) = prevX;
    }

    template<Scalar T>
    MDCell<T, 1> HamiltonMC<T>::makeDummyCell(VectorND<T> mass) {
        size_t dof = mass.getLength();
        return MDCell<T, 1>(dof, std::move(mass));
    }
}
