/*
 * Copyright 2023 Weibo He.
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

namespace Physica {
    template<Scalar T, unsigned int Dim = 3>
    class MonteCarlo {
        using MDCellType = MDCell<T, Dim>;
        using PositionMatrix = MDCellType::PositionMatrix;
        using VectorType = VectorND<T>;
        using MachineType = T::MachineType;

        MDCellType cell;
        PositionMatrix buffer;
        T repTemperature;
        T lastEnergy;
    public:
        MonteCarlo(MDCellType cell_, T temperatureT, T sigma);
        MonteCarlo(const MonteCarlo&) = default;
        MonteCarlo(MonteCarlo&&) noexcept = default;
        ~MonteCarlo() = default;
        /* Operators */
        MonteCarlo& operator=(MonteCarlo mc) noexcept;
        /* Operations */
        template<RNG R,
                 class ForceModel,
                 ExecutePolicy P>
        void nvt_step(const ForceModel& forceModel);
        void normalizePos() { cell.normalize(); }
        void swap(MonteCarlo& __restrict mc) noexcept;
        /* Getters */
        [[nodiscard]] const MDCellType& getCell() const noexcept { return cell; }
    };

    template<Scalar T, unsigned int Dim>
    MonteCarlo<T, Dim>::MonteCarlo(MDCellType cell_, T temperatureT, T sigma)
            : cell(std::move(cell_)), repTemperature(reciprocal(temperatureT)), lastEnergy(0) {
        buffer.resize(cell.getNumParticle(), Dim);
    }

    template<Scalar T, unsigned int Dim>
    MonteCarlo<T, Dim>& MonteCarlo<T, Dim>::operator=(MonteCarlo mc) noexcept {
        swap(mc);
        return *this;
    }

    template<Scalar T, unsigned int Dim>
    template<RNG R, class ForceModel, ExecutePolicy P>
    void MonteCarlo<T, Dim>::nvt_step(const ForceModel& forceModel) {
        for (size_t i = 0; i < cell.getNumParticle(); ++i) {
            for (unsigned int j = 0; j < Dim; ++j)
                buffer(i, j) = cell.getPos()(i, j) + T::template random_uniform<R>();
        }
        cell.swapPos(buffer);

        const T energy = forceModel.potentialV(cell);
        
        if (energy < lastEnergy) {
            lastEnergy = energy;
        }
        else {
            const T prob = exp((lastEnergy - energy) * repTemperature);
            const T temp = T::template random_uniform<R>();
            if (temp < prob)
                lastEnergy = energy;
            else
                cell.swapPos(buffer);
        }
    }

    template<Scalar T, unsigned int Dim>
    void MonteCarlo<T, Dim>::swap(MonteCarlo& __restrict mc) noexcept {
        assert(this != &mc && "[Error]: Self swap is likely a bug");
        cell.swap(mc.cell);
        buffer.swap(mc.buffer);
        repTemperature.swap(mc.repTemperature);
        lastEnergy.swap(mc.lastEnergy);
    }
}
