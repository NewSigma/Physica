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
    template<class ScalarType>
    class QuadraticSearch {
        using Vector3D = Vector3D<ScalarType>;

        Vector3D x;
        Vector3D y;
    public:
        QuadraticSearch(Vector3D x_, Vector3D y_);
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
        [[nodiscard]] ScalarType getOptimizedX() const noexcept { return x[getOptimalIndex()]; }
        [[nodiscard]] ScalarType getOptimizedY() const noexcept { return y[getOptimalIndex()]; }
    };

    template<class ScalarType>
    QuadraticSearch<ScalarType>::QuadraticSearch(Vector3D x_, Vector3D y_) : x(std::move(x_)), y(std::move(y_)) {}

    template<class ScalarType>
    template<class Functor>
    void QuadraticSearch<ScalarType>::step(Functor func) {
        using ResultType = typename std::invoke_result<Functor, ScalarType>::type;
        static_assert(std::is_same<ScalarType, ResultType>::value, "[Error]: Invalid functor");
        const ScalarType term1 = x[0] * (y[1] - y[2]);
        const ScalarType term2 = x[1] * (y[2] - y[0]);
        const ScalarType term3 = x[2] * (y[0] - y[1]);
        const ScalarType temp1 = x[0] * term1 + x[1] * term2 + x[2] * term3;
        const ScalarType temp2 = term1 + term2 + term3;
        const ScalarType new_x = temp1 / (temp2 * ScalarType(2));
        const ScalarType new_y = func(new_x);

        const int index = getOptimalIndex();
        x[index] = new_x;
        y[index] = new_y;
    }

    template<class ScalarType>
    void QuadraticSearch<ScalarType>::swap(QuadraticSearch& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        x.swap(obj.x);
        y.swap(obj.y);
    }

    template<class ScalarType>
    int QuadraticSearch<ScalarType>::getOptimalIndex() const noexcept {
        int index = 0;
        ScalarType max = std::numeric_limits<ScalarType>::min();
        for (int i = 0; i < 3; ++i) {
            if (y[i] > max) {
                index = i;
                max = y[i];
            }
        }
        return index;
    }
}
