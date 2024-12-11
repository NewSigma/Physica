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
    template<Scalar T>
    class device_obj<Loss<T>> {
        static_assert(!is_device_obj<T>::value, "[Error]: Nested device_obj<> is not allowed");
        using host_obj = Loss<T>;
    public:
        constexpr static bool IsTrainMode = T::isDiffable;
        using ValueType = T::ValueType;
        using LossType = std::conditional<IsTrainMode, device_obj<T>, T>::type;
    private:
        using ValueVector = VectorND<ValueType>;
        using DiffVector = Diff<ValueVector, DiffMode::Reverse, T::Order>;
        using VectorType = device_obj<typename std::conditional<IsTrainMode, DiffVector, ValueVector>::type>;
    public:
        [[nodiscard]] static LossType softmax(const VectorType& v, size_t label);
        [[nodiscard]] static LossType crossEntropy(const VectorType& v, size_t label);
    };

    template<Scalar T>
    device_obj<Loss<T>>::LossType
    device_obj<Loss<T>>::softmax(const VectorType& v, size_t label) {
        if constexpr (IsTrainMode) {
            const auto maximum = v.max();
            const auto temp = exp(v - maximum);
            return temp.calc(label) / temp.sum();
        }
        else
            return host_obj::softmax(v.toHost(), label);
    }

    template<Scalar T>
    device_obj<Loss<T>>::LossType
    device_obj<Loss<T>>::crossEntropy(const VectorType& v, size_t label) {
        assert(label < v.getLength() && "[Error]: The label is not exist");
        return -ln(softmax(v, label) + LossType(std::numeric_limits<ValueType>::min())); //Add minimum to avoid ln(0)
    }
}
