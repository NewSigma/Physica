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

#include <unordered_map>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

namespace Physica {
    template<Scalar T>
    class MomentSGD {
        static_assert(!Diffable<T>);

        using This = MomentSGD<T>;
        using VectorType = VectorND<T>;
        using MatrixType = DenseMatrix<T, MatrixOption::Col | MatrixOption::Element>;
    private:
        std::unordered_map<void*, std::variant<VectorType, MatrixType>> targetBufferMap;
        T lr;
        T moment;
    public:
        MomentSGD(T learnRate, T moment_);
        MomentSGD(const This&) = default;
        MomentSGD(This&&) noexcept = default;
        ~MomentSGD() = default;
        /* Operators */
        This& operator=(MomentSGD obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void step(Diffable auto& target);
        void step(Diffable auto& target, Diffable auto&... targets);
        void swap(This& __restrict obj) noexcept;
    };

    template<Scalar T>
    MomentSGD<T>::MomentSGD(T learnRate, T moment_)
            : lr(std::move(learnRate))
            , moment(std::move(moment_)) {
        assert(moment.isPositive() && "[Error]: Invalid moment");
    }

    template<Scalar T>
    void MomentSGD<T>::step(Diffable auto& target) {
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
    void MomentSGD<T>::step(Diffable auto& target, Diffable auto&... targets) {
        step(target);
        step(targets...);
    }

    template<Scalar T>
    void MomentSGD<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        lr.swap(obj.lr);
        targetBufferMap.swap(obj.targetBufferMap);
        moment.swap(obj.moment);
    }
}
