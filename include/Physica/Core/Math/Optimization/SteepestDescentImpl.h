/*
 * Copyright 2020-2023 WeiBo He.
 *
 * This file is part of Physica.

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
    template<class ScalarType, size_t Dim>
    SteepestDescent<ScalarType, Dim>::SteepestDescent(VectorType initial, ScalarType defaultStep_)
            : tryX(std::move(initial))
            , nowY(std::numeric_limits<ScalarType>::max())
            , defaultStep(defaultStep_)
            , stepSize(defaultStep_)
            , stepMultiple(0) {
        nowX = tryX;
    }

    template<class ScalarType, size_t Dim>
    SteepestDescent<ScalarType, Dim>& SteepestDescent<ScalarType, Dim>::operator=(SteepestDescent obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType, size_t Dim>
    template<class Functor>
    void SteepestDescent<ScalarType, Dim>::init(Functor func) {
        nowX = tryX;
        nowY = func(tryX);
        stepSize = defaultStep;
        stepMultiple = 0;
    }

    template<class ScalarType, size_t Dim>
    template<class Functor, class GradFunctor>
    bool SteepestDescent<ScalarType, Dim>::step(Functor func, GradFunctor grad) {
        VectorType normalGrad = grad(tryX);
        normalGrad.toUnit();
        tryX -= normalGrad * stepSize;
        ScalarType y = func(tryX);

        if (y < nowY) {
            nowX = tryX;
            nowY = y;
            stepSize *= 2.0;
            stepMultiple += 1;
            const bool stepIsTooLarge = stepMultiple == std::numeric_limits<unsigned char>::max();
            if (stepIsTooLarge) [[unlikely]]
                throw std::invalid_argument("Minimal cannot be achieved with a large step, are we dealing with a minimization problem?");
        }
        else {
            stepSize *= 0.5;
            stepMultiple -= 1;
            const bool isConverged = stepMultiple == 0;
            if (isConverged)
                return true;
        }
        return false;
    }

    template<class ScalarType, size_t Dim>
    template<class Functor, class GradFunctor>
    ScalarType SteepestDescent<ScalarType, Dim>::solve(Functor func, GradFunctor grad) {
        init(func);
        while (!step(func, grad));
        return nowY;
    }

    template<class ScalarType, size_t Dim>
    void SteepestDescent<ScalarType, Dim>::swap(SteepestDescent& obj) noexcept {
        tryX.swap(obj.tryX);
        nowX.swap(obj.nowX);
        nowY.swap(obj.nowY);
        defaultStep.swap(obj.defaultStep);
        stepSize.swap(obj.stepSize);
        std::swap(stepMultiple, obj.stepMultiple);
    }
}
