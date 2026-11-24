/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "../SteepestDescent.h"

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
    void SteepestDescent<T, Dim>::init(VectorType initial, std::invocable<VectorType> auto fn, std::invocable<VectorType> auto grad) {
        tryX = (std::move(initial));
        nowX = tryX;
        nowY = fn(tryX);
        gradG = grad(nowX);
    }

    template<Scalar T, size_t Dim>
    void SteepestDescent<T, Dim>::step(std::invocable<VectorType> auto fn, std::invocable<VectorType> auto grad) {
        VectorType direction = -gradG;
        const T stepSize = lineSearch.run(fn, grad, nowX, gradG, direction);

        nowX += direction * stepSize;
        nowY = fn(nowX);
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
    size_t SteepestDescent<T, Dim>::getDim() const noexcept {
        if constexpr (Dim == Dynamic)
            return gradG.getLength();
        return Dim;
    }
}
