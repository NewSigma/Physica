/*
 * Copyright 2023 WeiBo He.
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

#include "LayerBase.h"

namespace Physica::Core {
    template<class ScalarType> class Relu;

    namespace Internal {
        template<class T>
        class Traits<Relu<T>> {
        public:
            using ScalarType = T;
        }
    }

    template<class ScalarType>
    class Relu : public LayerBase<Relu<ScalarType>> {
        using Base = LayerBase<Relu<ScalarType>>;
        using VectorType = typename Base::VectorType;
    public:
        Relu() = default;
        Relu(const Relu&) = default;
        Relu(Relu&&) noexcept = default;
        ~Relu() = default;
        /* Operators */
        Relu& operator=(Relu obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] VectorType forward(const VectorType& x) const;
        inline void swap(Relu& obj) noexcept;
    };

    template<class ScalarType>
    typename Relu<ScalarType>::VectorType Relu<ScalarType>::forward(const VectorType& x) const {
        const size_t length = x.getLength();
        VectorType result(length);
        for (size_t i = 0; i < length; ++i)
            result[i] = x[i].isPositive() ? x[i] : ScalarType(0);
        return result;
    }

    template<class ScalarType>
    inline void Relu<ScalarType>::swap(Relu& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
    }
}
