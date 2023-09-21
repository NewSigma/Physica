/*
 * Copyright 2023 WeiBo He.
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

namespace Physica::Core {
    template<class ScalarType, unsigned int Dim = 3>
    class MonteCarlo {
        using MDCellType = MDCell<ScalarType, ScalarType, Dim>;
        using PositionMatrix = typename MDCellType::PositionMatrix;
        using VectorType = Vector<ScalarType>;
        using TrivialType = typename ScalarType::TrivialType;

        MDCellType cell;
        PositionMatrix buffer;
        ScalarType repTemperature;
        ScalarType lastEnergy;
        std::uniform_real_distribution<TrivialType> dist;
    public:
        MonteCarlo(MDCellType cell_, ScalarType temperatureT, ScalarType sigma);
        MonteCarlo(const MonteCarlo&) = default;
        MonteCarlo(MonteCarlo&&) noexcept = default;
        ~MonteCarlo() = default;
        /* Operators */
        MonteCarlo& operator=(MonteCarlo mc) noexcept;
        /* Operations */
        template<class RandomGenerator,
                 class ForceModel,
                 class Executor>
        void nvt_step(RandomGenerator& gen, const ForceModel& forceModel);
        void normalizePos() { cell.normalize(); }
        void swap(MonteCarlo& mc) noexcept;
        /* Getters */
        [[nodiscard]] const MDCellType& getCell() const noexcept { return cell; }
    };

    template<class ScalarType, unsigned int Dim>
    MonteCarlo<ScalarType, Dim>::MonteCarlo(MDCellType cell_, ScalarType temperatureT, ScalarType sigma)
            : cell(std::move(cell_)), repTemperature(reciprocal(temperatureT)), lastEnergy(0), dist(-TrivialType(sigma), TrivialType(sigma)) {
        buffer.resize(cell.getNumParticle(), Dim);
    }

    template<class ScalarType, unsigned int Dim>
    MonteCarlo<ScalarType, Dim>& MonteCarlo<ScalarType, Dim>::operator=(MonteCarlo mc) noexcept {
        swap(mc);
        return *this;
    }

    template<class ScalarType, unsigned int Dim>
    template<class RandomGenerator, class ForceModel, class Executor>
    void MonteCarlo<ScalarType, Dim>::nvt_step(RandomGenerator& gen, const ForceModel& forceModel) {
        std::uniform_real_distribution<> uniform_dist{};

        for (size_t i = 0; i < cell.getNumParticle(); ++i) {
            for (unsigned int j = 0; j < Dim; ++j)
                buffer(i, j) = cell.getPos()(i, j) + dist(gen);
        }
        cell.swapPos(buffer);

        const ScalarType energy = forceModel.potentialEnergy(cell);
        
        if (energy < lastEnergy) {
            lastEnergy = energy;
        }
        else {
            const ScalarType prob = exp((lastEnergy - energy) * repTemperature);
            const ScalarType temp = uniform_dist(gen);
            if (temp < prob)
                lastEnergy = energy;
            else
                cell.swapPos(buffer);
        }
    }

    template<class ScalarType, unsigned int Dim>
    void MonteCarlo<ScalarType, Dim>::swap(MonteCarlo& mc) noexcept {
        assert(this != &mc && "[Error]: Self swap is likely a bug");
        cell.swap(mc.cell);
        buffer.swap(mc.buffer);
        repTemperature.swap(mc.repTemperature);
        lastEnergy.swap(mc.lastEnergy);
        dist.swap(mc.dist);
    }
}
