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

namespace Physica {
    /**
     * Stochastic gradient descent for auto diff
     */
    template<Scalar T>
    class SGD {
        static_assert(!Diffable<T>);

        T learnRate;
    public:
        inline SGD(T lr);
        SGD(const SGD&) = default;
        SGD(SGD&&) noexcept = default;
        ~SGD() = default;
        /* Operators */
        SGD& operator=(SGD obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<Diffable U>
        void step(U& target) const;
        inline void swap(SGD& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] T getLearnRate() const noexcept { return learnRate; }
        /* Setters */
        inline void setLearnRate(T lr);
    };

    template<Scalar T>
    inline SGD<T>::SGD(T lr) {
        setLearnRate(std::move(lr));
    }

    template<Scalar T>
    template<Diffable U>
    void SGD<T>::step(U& target) const {
        if constexpr (Scalar<T>)
            target.value() -= learnRate * target.grad().value();
        else if constexpr (Vector<T> || Matrix<T>)
            target.values() -= learnRate * target.grads().values();
        else
            target.step(*this);
    }

    template<Scalar T>
    inline void SGD<T>::swap(SGD& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        learnRate.swap(obj.learnRate);
    }

    template<Scalar T>
    inline void SGD<T>::setLearnRate(T lr) {
        assert(lr.isPositive() && "[Error]: Invalid learn rate");
        learnRate = lr;
    }
}
