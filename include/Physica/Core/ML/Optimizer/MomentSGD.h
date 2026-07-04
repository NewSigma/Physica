/*
 * Copyright 2023-2026 Weibo He.
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

#include <any>
#include <unordered_map>
#include "OptBase.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

namespace Physica {
    template<Scalar T>
    class MomentSGD : public OptBase<T> {
        using This = MomentSGD<T>;
        using Base = OptBase<T>;

        template<class Value>
        struct Node;
    private:
        std::unordered_map<void*, std::any> targetNodeMap;
        T lr;
        T moment;
    public:
        MomentSGD(T learnRate, T moment_);
        MomentSGD(const This&) = default;
        MomentSGD(This&&) noexcept = default;
        ~MomentSGD() = default;
        /* Operators */
        This& operator=(This obj) noexcept;
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
    auto MomentSGD<T>::operator=(This obj) noexcept -> This& {
        swap(obj);
        return *this;
    }

    template<Scalar T>
    void MomentSGD<T>::step(Diffable auto& target) {
        using Target = decltype(target);
        using NodeT = Node<typename Base::template ValueT<Target>>;

        void* pTarget = (void*)(&target);
        const bool exist = targetNodeMap.contains(pTarget);
        auto& any = targetNodeMap[pTarget];
        if (!exist)
            any = NodeT(target);

        auto& lastGrad = std::any_cast<NodeT>(any).lastGrad;
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
        targetNodeMap.swap(obj.targetNodeMap);
        moment.swap(obj.moment);
    }

    template<Scalar T>
    template<class Value>
    struct MomentSGD<T>::Node {
        Value lastGrad;

        Node() = default;
        Node(const auto& src);
    };

    template<Scalar T>
    template<class Value>
    MomentSGD<T>::Node<Value>::Node(const auto& src) {
        if constexpr (Scalar<Value>)
            lastGrad = 0;
        else {
            lastGrad.resize(src);
            lastGrad.zeros();
        }
    }
}
