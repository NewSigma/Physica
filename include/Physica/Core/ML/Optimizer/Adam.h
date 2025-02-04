/*
 * Copyright 2021-2025 Weibo He.
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

#include <unordered_map>
#include "AdamImpl/AdamImpl.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] arXiv:1412.6980; https://doi.org/10.48550/arXiv.1412.6980
     */
    template<Scalar T>
    class Adam {
        static_assert(!Diffable<T>);
        using Args = AdamBase<T>::Args;

        Args args;
        std::unordered_map<void*, AdamBase<T>*> targetBufferMap;
    public:
        Adam(T lr = 1E-3, T beta1 = 0.9, T beta2 = 0.999, T epsilon = 1E-8, T decay = 0);
        Adam(const Adam&) = default;
        Adam(Adam&&) noexcept = default;
        ~Adam() = default;
        /* Operators */
        Adam& operator=(Adam obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class U>
        void step(U& target);

        void swap(Adam& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] Args& getArgs() const noexcept { return args; }
    };

    template<Scalar T>
    Adam<T>::Adam(T lr, T beta1, T beta2, T epsilon, T decay) {
        assert(lr.isPositive());
        assert(beta1.isPositive() && beta1 < T(1));
        assert(beta2.isPositive() && beta2 < T(1));
        assert(epsilon.isPositive());
        assert(!decay.isNegative() && decay < T(1));
        args.lr = lr;
        args.beta1 = beta1;
        args.beta2 = beta2;
        args.epsilon = epsilon;
        args.decay = decay;
    }

    template<Scalar T>
    template<class U>
    void Adam<T>::step(U& target) {
        void* pTarget = (void*)(&target);
        const bool exist = targetBufferMap.count(pTarget) != 0;
        auto& opt = targetBufferMap[pTarget];
        if (!exist) {
            if constexpr (Vector<U>)
                opt = new AdamImpl<U>(args, target.getLength());
            else {
                static_assert(Matrix<U>, "[Error]: Unexpected type");
                opt = new AdamImpl<U>(args, target.getRow(), target.getCol());
            }
        }
        opt->step(args, pTarget);
    }

    template<Scalar T>
    void Adam<T>::swap(Adam& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(args, obj.args);
        targetBufferMap.swap(obj.targetBufferMap);
    }
}
