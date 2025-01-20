/*
 * Copyright 2023-2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"
#include "Physica/Core/Scalar/Real.h"

namespace Physica::Core {
    /**
     * Stochastic gradient descent for auto diff
     */
    class SGD {
        float64 learnRate;
    public:
        inline SGD(float64 lr);
        SGD(const SGD&) = default;
        SGD(SGD&&) noexcept = default;
        ~SGD() = default;
        /* Operators */
        SGD& operator=(SGD obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<Diffable T>
        void step(T& obj) const;
        inline void swap(SGD& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] float64 getLearnRate() const noexcept { return learnRate; }
        /* Setters */
        inline void setLearnRate(float64 lr);
    };

    inline SGD::SGD(float64 lr) {
        setLearnRate(std::move(lr));
    }

    template<Diffable T>
    void SGD::step(T& obj) const {
        if constexpr (Scalar<T>)
            obj.value() -= learnRate * obj.grad().value();
        else if constexpr (Vector<T> || Matrix<T>)
            obj.values() -= learnRate * obj.grads().values();
        else
            obj.step(*this);
    }

    inline void SGD::swap(SGD& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        learnRate.swap(obj.learnRate);
    }

    inline void SGD::setLearnRate(float64 lr) {
        assert(!lr.isPositive() && "[Error]: 0 learn rate does nothing");
        learnRate = lr;
    }
}
