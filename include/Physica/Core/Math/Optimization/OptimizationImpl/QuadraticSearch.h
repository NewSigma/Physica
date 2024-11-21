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

namespace Physica::Core {
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
        template<class Functor>
        void step(Functor func);
        void swap(QuadraticSearch& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] int getOptimalIndex() const noexcept;
        [[nodiscard]] T getOptimizedX() const noexcept { return x[getOptimalIndex()]; }
        [[nodiscard]] T getOptimizedY() const noexcept { return y[getOptimalIndex()]; }
    };

    template<Scalar T>
    QuadraticSearch<T>::QuadraticSearch(Vector3D<T> x_, Vector3D<T> y_) : x(std::move(x_)), y(std::move(y_)) {}

    template<Scalar T>
    template<class Functor>
    void QuadraticSearch<T>::step(Functor func) {
        using ResultType = std::invoke_result<Functor, T>::type;
        static_assert(std::is_same<T, ResultType>::value, "[Error]: Invalid functor");
        const T term1 = x[0] * (y[1] - y[2]);
        const T term2 = x[1] * (y[2] - y[0]);
        const T term3 = x[2] * (y[0] - y[1]);
        const T temp1 = x[0] * term1 + x[1] * term2 + x[2] * term3;
        const T temp2 = term1 + term2 + term3;
        const T new_x = temp1 / (temp2 * T(2));
        const T new_y = func(new_x);

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

    template<Scalar T>
    int QuadraticSearch<T>::getOptimalIndex() const noexcept {
        int index = 0;
        T max = std::numeric_limits<T>::min();
        for (int i = 0; i < 3; ++i) {
            if (y[i] > max) {
                index = i;
                max = y[i];
            }
        }
        return index;
    }
}
