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

#include "MomentSGD.h"

namespace Physica {
    template<Scalar T>
    class device_obj<MomentSGD<T>> {
        static_assert(!Diffable<T>);

        using This = device_obj<MomentSGD<T>>;
        using VectorType = device_obj<VectorND<T>>;
        using MatrixType = device_obj<DenseMatrix<T>>;
    private:
        std::unordered_map<void*, std::variant<VectorType, MatrixType>> targetBufferMap;
        T lr;
        T moment;
    public:
        device_obj(T learnRate, T moment_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void step(Diffable auto& target);
        void step(Diffable auto& target, Diffable auto&... targets);
        void swap(This& __restrict obj) noexcept;
    };

    template<Scalar T>
    device_obj<MomentSGD<T>>::device_obj(T learnRate, T moment_)
            : lr(std::move(learnRate))
            , moment(std::move(moment_)) {
        assert(moment.isPositive() && "[Error]: Invalid moment");
    }

    template<Scalar T>
    void device_obj<MomentSGD<T>>::step(Diffable auto& target) {
        using U = decltype(target);
        using BufferType = std::conditional<Vector<U>, VectorType, MatrixType>::type;

        void* pTarget = (void*)(&target);
        const bool exist = targetBufferMap.count(pTarget) != 0;
        auto& var = targetBufferMap[pTarget];
        if (!exist) {
            if constexpr (Vector<U>)
                var = BufferType(target.getLength());
            else {
                static_assert(Matrix<U>, "[Error]: Unexpected type");
                var = BufferType(target.getRow(), target.getCol());
            }
            std::get<BufferType>(var).zeros();
        }

        auto& lastGrad = std::get<BufferType>(var);
        lastGrad = moment * lastGrad + target.grads();
        target.values() -= lr * lastGrad;
    }

    template<Scalar T>
    void device_obj<MomentSGD<T>>::step(Diffable auto& target, Diffable auto&... targets) {
        step(target);
        step(targets...);
    }

    template<Scalar T>
    void device_obj<MomentSGD<T>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        targetBufferMap.swap(obj.targetBufferMap);
        lr.swap(obj.lr);
        moment.swap(obj.moment);
    }
}
