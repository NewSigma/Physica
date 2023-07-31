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

#include <limits>
#include "Scalar.h"

namespace Physica::Core {
    template<class ScalarType> class Differentiable;
    /**
     * Auto differential support for scalars
     */
    template<class ScalarType>
    class Differentiable {
        ScalarType value;
        ScalarType tangent;
    public:
        Differentiable(ScalarType value_);
        Differentiable(ScalarType value_, ScalarType tangent_);
        Differentiable(const Differentiable&) = default;
        Differentiable(Differentiable&&) noexcept = default;
        ~Differentiable() = default;
        /* Operators */
        Differentiable& operator=(Differentiable obj) noexcept;
        [[nodiscard]] inline Differentiable operator-() const;
        /* Operations */
        void swap(Differentiable& obj) noexcept;
        /* Getters */
        [[nodiscard]] const ScalarType& getValue() const noexcept { return value; }
        [[nodiscard]] const ScalarType& getTangent() const noexcept { return tangent; }
    };

    template<class ScalarType>
    Differentiable<ScalarType>::Differentiable(ScalarType value_)
        : value(std::move(value_)), tangent(0) {}

    template<class ScalarType>
    Differentiable<ScalarType>::Differentiable(ScalarType value_, ScalarType tangent_)
        : value(std::move(value_)), tangent(std::move(tangent_)) {}

    template<class ScalarType>
    Differentiable<ScalarType>& Differentiable<ScalarType>::operator=(Differentiable obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType>
    inline Differentiable<ScalarType> Differentiable<ScalarType>::operator-() const {
        return {-value, -tangent};
    }

    template<class ScalarType>
    void Differentiable<ScalarType>::swap(Differentiable& obj) noexcept {
        value.swap(obj.value);
        tangent.swap(obj.tangent);
    }
}

namespace std {
    template<class ScalarType>
    struct numeric_limits<Physica::Core::Differentiable<ScalarType>> : public numeric_limits<ScalarType> {};
}

#include "DifferentiableImpl/ElementaryFunction.h"
