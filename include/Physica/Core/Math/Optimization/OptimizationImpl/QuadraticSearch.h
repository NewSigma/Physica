/*
 * Copyright 2023-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica {
    template<Scalar T>
    class QuadraticSearch {
        Vector3D<T> x;
        Vector3D<T> y;
    public:
        QuadraticSearch(Vector3D<T> x_, Vector3D<T> y_);
        QuadraticSearch(const QuadraticSearch&) = default;
        QuadraticSearch(QuadraticSearch&&) noexcept = default;
        ~QuadraticSearch() = default;
        /* Operators */
        QuadraticSearch& operator=(QuadraticSearch obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void step(std::invocable<T> auto fn);
        void swap(QuadraticSearch& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getOptimalIndex() const noexcept { return y.argmax(); }
        [[nodiscard]] T getOptimizedX() const noexcept { return x[getOptimalIndex()]; }
        [[nodiscard]] T getOptimizedY() const noexcept { return y[getOptimalIndex()]; }
    };

    template<Scalar T>
    QuadraticSearch<T>::QuadraticSearch(Vector3D<T> x_, Vector3D<T> y_) : x(std::move(x_)), y(std::move(y_)) {}

    template<Scalar T>
    void QuadraticSearch<T>::step(std::invocable<T> auto fn) {
        const T term1 = x[0] * (y[1] - y[2]);
        const T term2 = x[1] * (y[2] - y[0]);
        const T term3 = x[2] * (y[0] - y[1]);
        const T temp1 = x[0] * term1 + x[1] * term2 + x[2] * term3;
        const T temp2 = term1 + term2 + term3;
        const T new_x = temp1 / (temp2 * T(2));
        const T new_y = fn(new_x);

        const int index = getOptimalIndex();
        x[index] = new_x;
        y[index] = new_y;
    }

    template<Scalar T>
    void QuadraticSearch<T>::swap(QuadraticSearch& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        x.swap(obj.x);
        y.swap(obj.y);
    }
}
