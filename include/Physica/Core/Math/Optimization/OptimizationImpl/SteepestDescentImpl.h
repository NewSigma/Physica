/*
 * Copyright 2020-2023 Weibo He.
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

namespace Physica {
    template<Scalar T, size_t Dim>
    SteepestDescent<T, Dim>::SteepestDescent(T maxStepSize, T decreaseCondNum, T curvatureCondNum)
            : nowY(std::numeric_limits<T>::max())
            , lineSearch(maxStepSize, decreaseCondNum, curvatureCondNum) {}

    template<Scalar T, size_t Dim>
    SteepestDescent<T, Dim>& SteepestDescent<T, Dim>::operator=(SteepestDescent obj) noexcept {
        swap(obj);
        return *this;
    }

    template<Scalar T, size_t Dim>
    template<class Functor, class GradFunctor>
    void SteepestDescent<T, Dim>::init(VectorType initial, Functor func, GradFunctor grad) {
        tryX = (std::move(initial));
        nowX = tryX;
        nowY = func(tryX);
        gradG = grad(nowX);
    }

    template<Scalar T, size_t Dim>
    template<class Functor, class GradFunctor>
    void SteepestDescent<T, Dim>::step(Functor func, GradFunctor grad) {
        VectorType direction = -gradG;
        const T stepSize = lineSearch.run(func, grad, nowX, gradG, direction);

        nowX += direction * stepSize;
        nowY = func(nowX);
        gradG = grad(nowX);
    }

    template<Scalar T, size_t Dim>
    void SteepestDescent<T, Dim>::swap(SteepestDescent& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        gradG.swap(obj.gradG);
        tryX.swap(obj.tryX);
        nowX.swap(obj.nowX);
        nowY.swap(obj.nowY);
        lineSearch.swap(obj.lineSearch);
    }

    template<Scalar T, size_t Dim>
    inline size_t SteepestDescent<T, Dim>::getDim() const noexcept {
        if constexpr (Dim == Dynamic)
            return gradG.getLength();
        return Dim;
    }
}
