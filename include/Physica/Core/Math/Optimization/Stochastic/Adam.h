/*
 * Copyright 2021-2024 Weibo He.
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

namespace Physica::Core {
    /**
     * Reference:
     * [1] Adam: A Method for Stochastic Optimization (arXiv:1412.6980 [cs.LG])
     */
    template<Scalar T, Vector U>
    class Adam {
        const Array<T, 6> args;
        U params;
    public:
        Adam(const Array<T, 6>& args_);
        ~Adam() = default;
        /* Operations */
        template<class Function>
        void compute(Function func, const U& params_, size_t maxIteration);
        /* Getters */
        [[nodiscard]] const U& getParams() const noexcept { return params; }
    private:
        template<class Function>
        [[nodiscard]] U gradient(Function func);
    };
    /**
     * Five args:
     * args[0]: stepSize, real number in (0, inf)
     * args[1]: real number in [0, 1)
     * args[2]: real number in [0, 1)
     * args[3]: Exponential decay rates for the moment estimates, real number in [0, 1)
     * args[4]: Very small real number in (0, inf)
     * args[5]: Expected relative error, small real number in (0, inf)
     */
    template<Scalar T, Vector U>
    Adam<T, U>::Adam(const Array<T, 6>& args_) : args(args_), params() {
        assert(args[0].isPositive());
        assert(args[1].isPositive() && args[1] < T(1));
        assert(args[2].isPositive() && args[2] < T(1));
        assert(args[3].isPositive() && args[3] < T(1));
        assert(args[4].isPositive());
        assert(args[5].isPositive());
    }
    /**
     * \class Function is declared like:
     * T func(const U& params)
     * 
     * \param maxIteration
     * Set to 0 to disable this criteria
     */
    template<Scalar T, Vector U>
    template<class Function>
    void Adam<T, U>::compute(Function func, const U& params_, size_t maxIteration) {
        params = params_;
        auto m = U::zeros(params.getLength());
        auto v = U::zeros(params.getLength());
        T beta1 = args[1];
        size_t count = 0;
        U temp(params.getLength());
        bool stop = false;
        do {
            U g = gradient(func);
            const T beta1_1 = T(1) - beta1;
            m = m * beta1 + g * beta1_1;
            v = v * args[2] + hadamard(g, g) * T(T(1) -  args[2]);
            const T alpha = args[0] / beta1_1 * sqrt(T(T(1) - args[2]));
            temp = params - alpha * hadamard(m, reciprocal(sqrt(v) + args[4]));
            beta1 *= args[3];
            ++count;
            const bool meetRelativeCriteria = vectorNear(params, temp, double(args[5]));
            stop = (maxIteration != 0 && count > maxIteration) || meetRelativeCriteria;
            params = temp;
        } while(!stop);
    }

    template<Scalar T, Vector U>
    template<class Function>
    U Adam<T, U>::gradient(Function func) {
        const size_t length = params.getLength();
        U result(length);
        const T y0 = func(std::cref(params));
        for (size_t i = 0; i < length; ++i) {
            T new_param = params[i] + args[0];
            params[i].swap(new_param);
            const T y1 = func(std::cref(params));
            result[i] = (y1 - y0) / args[0];
            params[i].swap(new_param);
        }
        return result;
    }
}
