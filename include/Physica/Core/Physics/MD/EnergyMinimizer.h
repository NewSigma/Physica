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

#include "MDCell.h"
#include "Physica/Core/Math/Calculus/Differential.h"
#include "Physica/Core/Math/Optimization/SteepestDescent.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/ReshapedVector.h"

namespace Physica::Core {
    template<class ScalarType, unsigned int Dim = 3>
    class EnergyMinimizer {
    public:
        using MDCellType = MDCell<ScalarType, Dim>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using VectorType = Vector<ScalarType, Dynamic>;
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
        template<class ForceModel, class Executor, class Optimizer>
        void pre_pos_step(const ForceModel& model, Optimizer& optimizer);
        template<class ForceModel, class Executor, class Optimizer>
        void pos_step(const ForceModel& model, Optimizer& optimizer);
        template<class ForceModel, class Executor, class Optimizer>
        void pre_lattice2D_step(const ForceModel& model, Optimizer& optimizer, ScalarType diffStep);
        template<class ForceModel, class Executor, class Optimizer>
        void lattice2D_step(const ForceModel& model, Optimizer& optimizer, ScalarType diffStep);
        void swap(EnergyMinimizer& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const MDCellType& getCell() const noexcept { return cell; }
    private:
        template<class ForceModel> auto makePosStepFunc(const ForceModel& model);
        template<class ForceModel, class Executor> auto makePosStepGrad(const ForceModel& model);
        template<class ForceModel> auto makeLattice2DStepFunc(const ForceModel& model);
        template<class ForceModel> auto makeLattice2DStepGrad(const ForceModel& model, ScalarType diffStep);
    };

    template<class ScalarType, unsigned int Dim>
    EnergyMinimizer<ScalarType, Dim>::EnergyMinimizer(MDCellType cell_)
            : cell(std::move(cell_)) {}

    template<class ScalarType, unsigned int Dim>
    EnergyMinimizer<ScalarType, Dim>&
    EnergyMinimizer<ScalarType, Dim>::operator=(EnergyMinimizer<ScalarType, Dim> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, unsigned int Dim>
    template<class ForceModel, class Executor, class Optimizer>
    void EnergyMinimizer<ScalarType, Dim>::pre_pos_step(const ForceModel& model, Optimizer& optimizer) {
        const auto func = makePosStepFunc(model);
        const auto grad = makePosStepGrad<ForceModel, Executor>(model);
        optimizer.init(cell.getPos().flatten(), func, grad);
    }

    template<class ScalarType, unsigned int Dim>
    template<class ForceModel, class Executor, class Optimizer>
    void EnergyMinimizer<ScalarType, Dim>::pos_step(const ForceModel& model, Optimizer& optimizer) {
        assert(optimizer.getDim() == cell.getDOF() && "[Error]: pre_pos_step must be called before this function");
        const auto func = makePosStepFunc(model);
        const auto grad = makePosStepGrad<ForceModel, Executor>(model);
        optimizer.step(func, grad);
        cell.setPos(optimizer.getArgX().reshape(cell.getPos()));
    }

    template<class ScalarType, unsigned int Dim>
    template<class ForceModel, class Executor, class Optimizer>
    void EnergyMinimizer<ScalarType, Dim>::pre_lattice2D_step(const ForceModel& model, Optimizer& optimizer, ScalarType diffStep) {
        assert(cell.getLattice()(0, 1).isZero());
        assert(cell.getLattice()(0, 2).isZero());
        assert(cell.getLattice()(1, 2).isZero());
        assert(cell.getLattice()(2, 0).isZero());
        assert(cell.getLattice()(2, 1).isZero());
        const auto func = makeLattice2DStepFunc(model);
        const auto grad = makeLattice2DStepGrad(model, diffStep);
        const auto& lattice = cell.getLattice();
        optimizer.init({lattice(0, 0), lattice(1, 0), lattice(1, 1)}, func, grad);
    }

    template<class ScalarType, unsigned int Dim>
    template<class ForceModel, class Executor, class Optimizer>
    void EnergyMinimizer<ScalarType, Dim>::lattice2D_step(const ForceModel& model, Optimizer& optimizer, ScalarType diffStep) {
        assert(optimizer.getDim() == 3 && "[Error]: pre_lattice2D_step must be called before this function");
        const auto func = makeLattice2DStepFunc(model);
        const auto grad = makeLattice2DStepGrad(model, diffStep);
        optimizer.step(func, grad);

        LatticeMatrix lattice = cell.getLattice();
        const auto& argX = optimizer.getArgX();
        lattice(0, 0) = argX[0];
        lattice(1, 0) = argX[1];
        lattice(1, 1) = argX[2];
        cell.setLattice(lattice);
    }

    template<class ScalarType, unsigned int Dim>
    void EnergyMinimizer<ScalarType, Dim>::swap(EnergyMinimizer& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        cell.swap(obj.cell);
    }

    template<class ScalarType, unsigned int Dim>
    template<class ForceModel>
    auto EnergyMinimizer<ScalarType, Dim>::makePosStepFunc(const ForceModel& model) {
        auto func = [this, &model](const VectorType& v) -> ScalarType {
            const auto temp = MDCellType(cell.getLattice(), v.reshape(cell.getPos()), cell.getMassVec());
            return model.potentialV(temp);
        };
        return func;
    }

    template<class ScalarType, unsigned int Dim>
    template<class ForceModel, class Executor>
    auto EnergyMinimizer<ScalarType, Dim>::makePosStepGrad(const ForceModel& model) {
        const auto grad = [this, &model](const VectorType& v) -> VectorType {
            const auto temp = MDCellType(cell.getLattice(), v.reshape(cell.getPos()), cell.getMassVec());
            return -model.template force<Executor, true>(temp);
        };
        return grad;
    }

    template<class ScalarType, unsigned int Dim>
    template<class ForceModel>
    auto EnergyMinimizer<ScalarType, Dim>::makeLattice2DStepFunc(const ForceModel& model) {
        const auto func = [this, &model](const VectorType& v) -> ScalarType {
            LatticeMatrix lattice = cell.getLattice();
            lattice(0, 0) = v[0];
            lattice(1, 0) = v[1];
            lattice(1, 1) = v[2];
            const auto temp = MDCellType(std::move(lattice), cell.getPos(), cell.getMassVec());
            return model.potentialV(temp);
        };
        return func;
    }

    template<class ScalarType, unsigned int Dim>
    template<class ForceModel>
    auto EnergyMinimizer<ScalarType, Dim>::makeLattice2DStepGrad(const ForceModel& model, ScalarType diffStep) {
        const auto func = makeLattice2DStepFunc(model);
        const auto grad = [func, diffStep](const VectorType& v) -> VectorType {
            VectorType result(3);
            for (int i = 0; i < int(result.getLength()); ++i) {
                result[i] = Differential<ScalarType>::doublePoint([i, func, &v](ScalarType x) {
                    VectorType v1 = v;
                    v1[i] = x;
                    return func(v1);
                }, v[i], diffStep);
            }
            return result;
        };
        return grad;
    }
}
