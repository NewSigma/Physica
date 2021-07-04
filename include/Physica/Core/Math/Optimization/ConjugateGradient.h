/*
 * Copyright 2021 WeiBo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector.h"

namespace Physica::Core::Math {
    template<class Scalar, class Function, class Vector>
    class ConjugateGradient {
        Function func;
        Vector x;
        Scalar epsilon;
        Scalar minStepSize;
        Scalar minimal;
    public:
        ConjugateGradient(Function func_, Vector x_, Scalar epsilon_, Scalar minStepSize_);
        ~ConjugateGradient() = default;
        /* Operations */
        Scalar compute();
    private:
        [[nodiscard]] Vector gradient();
        [[nodiscard]] Scalar climbHill(const Vector& direction);
    };

    template<class Scalar, class Function, class Vector>
    ConjugateGradient<Scalar, Function, Vector>::ConjugateGradient(Function func_,
                                                                    Vector x_,
                                                                    Scalar epsilon_,
                                                                    Scalar minStepSize_)
            : func(func_), x(std::move(x_)), epsilon(epsilon_), minStepSize(minStepSize_) {
        assert(minStepSize.isPositive());
    }

    template<class Scalar, class Function, class Vector>
    Scalar ConjugateGradient<Scalar, Function, Vector>::compute() {
        const auto length = x.getLength();

        Scalar y_temp = func(std::cref(x));
        Vector gradient_temp = gradient();
        Vector search_direction = -gradient_temp;
        for (size_t i = 0; i <= length; ++i) {
            const Scalar factor = climbHill(search_direction);
            x += factor * search_direction;

            Scalar y_new = func(std::cref(x));
            const Scalar delta = abs(Scalar(y_temp - y_new));
            y_temp = std::move(y_new);

            if (delta < epsilon)
                break;

            Vector gradient_new = gradient();
            if (i == length) {
                i = 0;
                search_direction = -gradient_new;
                continue;
            }

            const Scalar norm = search_direction.norm();
            const Scalar alpha = gradient_new.squaredNorm() / gradient_temp.squaredNorm();
            search_direction = norm * alpha * search_direction - gradient_new;
            gradient_temp = std::move(gradient_new);

            if (!(gradient_temp * search_direction).isNegative()) {
                i = 0;
                search_direction = -gradient_temp;
                continue;
            }
        }
        return y_temp;
    }

    template<class Scalar, class Function, class Vector>
    Vector ConjugateGradient<Scalar, Function, Vector>::gradient() {
        const size_t length = x.getLength();
        Vector result(length);
        const Scalar y0 = func(std::cref(x));
        for (size_t i = 0; i < length; ++i) {
            Scalar new_x = x[i] + minStepSize;
            std::swap(x[i], new_x);
            const Scalar y1 = func(std::cref(x));
            result[i] = (y1 - y0) / minStepSize;
            std::swap(x[i], new_x);
        }
        return result;
    }

    template<class Scalar, class Function, class Vector>
    Scalar ConjugateGradient<Scalar, Function, Vector>::climbHill(const Vector& direction) {
        Scalar y0 = func(std::cref(x));
        Scalar y;
        Scalar result = minStepSize;
        bool stop;
        do {
            y = func(x + result * direction);
            stop = y > y0;
            result += minStepSize;
            y0 = std::move(y);
        } while(!stop);
        return result;
    }
}
