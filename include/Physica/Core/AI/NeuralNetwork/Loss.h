/*
 * Copyright 2023-2024 Weibo He.
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
    template<Scalar T>
    class Loss {
        static_assert(!is_device_obj<T>::value, "[Error]: Include corresponding *.cuh file to enable CUDA support");
    public:
        constexpr static bool IsTrainMode = T::isDiffable;
        using ValueType = T::ValueType;
        using LossType = T;
    private:
        using VectorType = VectorND<T>;
    public:
        [[nodiscard]] static T softmax(const VectorType& v, size_t label);
        [[nodiscard]] static T crossEntropy(const VectorType& v, size_t label);
        [[nodiscard]] static T polar_rate(T train_loss, T valid_loss);
        [[nodiscard]] static T mixed_loss(T train_loss, T valid_loss);
    };

    template<Scalar T>
    T Loss<T>::softmax(const VectorType& v, size_t label) {
        const T maximum = v.max();
        const T factor = reciprocal(exp(v - maximum).sum());
        return exp(v.calc(label) - maximum) * factor;
    }

    template<Scalar T>
    T Loss<T>::crossEntropy(const VectorType& v, size_t label) {
        assert(label < v.getLength() && "[Error]: The label is not exist");
        return -ln(softmax(v, label) + T(std::numeric_limits<ValueType>::min())); //Add minimum to avoid ln(0)
    }
    /**
     * \returns polarization rate, the lower the better, minus value means overfitting
     */
    template<Scalar T>
    T Loss<T>::polar_rate(T train_loss, T valid_loss) {
        const T total = train_loss + valid_loss;
        const T delta = train_loss - valid_loss;
        if (abs(delta) * T(std::numeric_limits<T>::epsilon()) >= total)
            return T(0);
        return delta / total * T(2);
    }

    template<Scalar T>
    T Loss<T>::mixed_loss(T train_loss, T valid_loss) {
        return std::max(train_loss, valid_loss) + abs(train_loss - valid_loss);
    }
}
