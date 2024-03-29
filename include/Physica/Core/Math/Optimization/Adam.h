/*
 * Copyright 2021-2022 WeiBo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"
#include "Physica/Utils/TestHelper.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] Adam: A Method for Stochastic Optimization (arXiv:1412.6980 [cs.LG])
     */
    template<class ScalarType, class VectorType>
    class Adam {
        const Utils::Array<ScalarType, 6> args;
        VectorType params;
    public:
        Adam(const Utils::Array<ScalarType, 6>& args_);
        ~Adam() = default;
        /* Operations */
        template<class Function>
        void compute(Function func, const VectorType& params_, size_t maxIteration);
        /* Getters */
        [[nodiscard]] const VectorType& getParams() const noexcept { return params; }
    private:
        template<class Function>
        [[nodiscard]] VectorType gradient(Function func);
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
    template<class ScalarType, class VectorType>
    Adam<ScalarType, VectorType>::Adam(const Utils::Array<ScalarType, 6>& args_) : args(args_), params() {
        assert(args[0].isPositive());
        assert(args[1].isPositive() && args[1] < ScalarType::One());
        assert(args[2].isPositive() && args[2] < ScalarType::One());
        assert(args[3].isPositive() && args[3] < ScalarType::One());
        assert(args[4].isPositive());
        assert(args[5].isPositive());
    }
    /**
     * \class Function is declared like:
     * ScalarType func(const VectorType& params)
     * 
     * \param maxIteration
     * Set to 0 to disable this criteria
     */
    template<class ScalarType, class VectorType>
    template<class Function>
    void Adam<ScalarType, VectorType>::compute(Function func, const VectorType& params_, size_t maxIteration) {
        params = params_;
        VectorType m = VectorType::Zeros(params.getLength());
        VectorType v = VectorType::Zeros(params.getLength());
        ScalarType beta1 = args[1];
        size_t count = 0;
        VectorType temp(params.getLength());
        bool stop = false;
        do {
            VectorType g = gradient(func);
            const ScalarType beta1_1 = ScalarType::One() - beta1;
            m = m * beta1 + g * beta1_1;
            v = v * args[2] + hadamard(g, g) * ScalarType(ScalarType::One() -  args[2]);
            const ScalarType alpha = args[0] / beta1_1 * sqrt(ScalarType(ScalarType::One() - args[2]));
            temp = params - alpha * hadamard(m, reciprocal(sqrt(v) + args[4]));
            beta1 *= args[3];
            ++count;
            const bool meetRelativeCriteria = Utils::vectorNear(params, temp, double(args[5]));
            stop = (maxIteration != 0 && count > maxIteration) || meetRelativeCriteria;
            params = temp;
        } while(!stop);
    }

    template<class ScalarType, class VectorType>
    template<class Function>
    VectorType Adam<ScalarType, VectorType>::gradient(Function func) {
        using std::swap;
        const size_t length = params.getLength();
        VectorType result(length);
        const ScalarType y0 = func(std::cref(params));
        for (size_t i = 0; i < length; ++i) {
            ScalarType new_param = params[i] + args[0];
            swap(params[i], new_param);
            const ScalarType y1 = func(std::cref(params));
            result[i] = (y1 - y0) / args[0];
            swap(params[i], new_param);
        }
        return result;
    }
}
