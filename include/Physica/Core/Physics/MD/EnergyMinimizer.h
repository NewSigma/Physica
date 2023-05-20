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

#include "MDCell.h"
#include "Physica/Core/Math/Optimization/SteepestDescent.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/ReshapedVector.h"

namespace Physica::Core {
    template<class ScalarType, class PosScalarType, unsigned int Dim = 3>
    class EnergyMinimizer {
    public:
        using MDCellType = MDCell<ScalarType, PosScalarType, Dim>;
        using VectorType = Vector<PosScalarType, Dynamic>;
    private:
        MDCellType cell;
    public:
        EnergyMinimizer(MDCellType cell_);
        EnergyMinimizer(const EnergyMinimizer&) = default;
        EnergyMinimizer(EnergyMinimizer&&) noexcept = default;
        ~EnergyMinimizer() = default;
        /* Operators */
        EnergyMinimizer& operator=(EnergyMinimizer obj) noexcept;
        /* Operations */
        template<class ForceModel, class Optimizer>
        void init(const ForceModel& model, Optimizer& optimizer);
        template<class ForceModel, class Executor, class Optimizer>
        bool pos_step(const ForceModel& model, Optimizer& optimizer);
        void swap(EnergyMinimizer& obj) noexcept;
        /* Getters */
        [[nodiscard]] const MDCellType& getCell() const noexcept { return cell; }
    };

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    EnergyMinimizer<ScalarType, PosScalarType, Dim>::EnergyMinimizer(MDCellType cell_)
            : cell(std::move(cell_)) {}

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    EnergyMinimizer<ScalarType, PosScalarType, Dim>&
    EnergyMinimizer<ScalarType, PosScalarType, Dim>::operator=(EnergyMinimizer<ScalarType, PosScalarType, Dim> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class ForceModel, class Optimizer>
    void EnergyMinimizer<ScalarType, PosScalarType, Dim>::init(const ForceModel& model, Optimizer& optimizer) {
        optimizer.init(cell.getPos().flatten(), [this, &model](const VectorType& v) -> ScalarType {
            const auto temp = MDCellType(cell.getLattice(), v.reshape(cell.getPos()), cell.getMassVec());
            return model.potentialEnergy(temp);
        });
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    template<class ForceModel, class Executor, class Optimizer>
    bool EnergyMinimizer<ScalarType, PosScalarType, Dim>::pos_step(const ForceModel& model, Optimizer& optimizer) {
        const auto func = [this, &model](const VectorType& v) -> ScalarType {
            const auto temp = MDCellType(cell.getLattice(), v.reshape(cell.getPos()), cell.getMassVec());
            return model.potentialEnergy(temp);
        };
        const auto grad = [this, &model](const VectorType& v) -> VectorType {
            const auto temp = MDCellType(cell.getLattice(), v.reshape(cell.getPos()), cell.getMassVec());
            return -model.template force<Executor, true>(temp);
        };
        const bool isDone = optimizer.step(func, grad);
        cell.setPos(optimizer.getArgX().reshape(cell.getPos()));
        return isDone;
    }

    template<class ScalarType, class PosScalarType, unsigned int Dim>
    void EnergyMinimizer<ScalarType, PosScalarType, Dim>::swap(EnergyMinimizer& obj) noexcept {
        cell.swap(obj.cell);
    }
}
