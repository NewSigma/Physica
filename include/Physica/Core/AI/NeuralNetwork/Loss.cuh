/*
 * Copyright 2024 Weibo He.
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

#include "Loss.h"

namespace Physica::Core {
    template<class ScalarType>
    class device_obj<Loss<ScalarType>> {
        static_assert(!is_device_obj<ScalarType>::value, "[Error]: Nested device_obj<> is not allowed");
        using host_obj = Loss<ScalarType>;
    public:
        constexpr static bool IsTrainMode = ScalarType::isDifferentiable;
        using PlainScalar = typename ScalarType::PlainScalar;
        using LossType = typename std::conditional<IsTrainMode, device_obj<ScalarType>, ScalarType>::type;
    private:
        using PlainVector = Vector<PlainScalar>;
        using DiffVector = Diff<PlainVector, DiffMode::Reverse, ScalarType::Order>;
        using VectorType = device_obj<typename std::conditional<IsTrainMode, DiffVector, PlainVector>::type>;
    public:
        [[nodiscard]] static LossType softmax(const VectorType& v, size_t label);
        [[nodiscard]] static LossType crossEntropy(const VectorType& v, size_t label);
    };

    template<class ScalarType>
    typename device_obj<Loss<ScalarType>>::LossType
    device_obj<Loss<ScalarType>>::softmax(const VectorType& v, size_t label) {
        if constexpr (IsTrainMode) {
            const auto maximum = v.max();
            const auto temp = exp(v - maximum);
            return temp.calc(label) / temp.sum();
        }
        else
            return host_obj::softmax(v.toHost(), label);
    }

    template<class ScalarType>
    typename device_obj<Loss<ScalarType>>::LossType
    device_obj<Loss<ScalarType>>::crossEntropy(const VectorType& v, size_t label) {
        assert(label < v.getLength() && "[Error]: The label is not exist");
        return -ln(softmax(v, label) + LossType(std::numeric_limits<PlainScalar>::min())); //Add minimum to avoid ln(0)
    }
}
