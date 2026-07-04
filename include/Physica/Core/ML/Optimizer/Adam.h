/*
 * Copyright 2021-2026 Weibo He.
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
    /**
     * Reference:
     * [1] arXiv:1412.6980; https://doi.org/10.48550/arXiv.1412.6980
     * [2] pytorch; https://pytorch.org/docs/stable/generated/torch.optim.Adam.html
     */
    template<Scalar T>
    class Adam : public OptBase<T> {
        using This = Adam<T>;
        using Base = OptBase<T>;

        template<class Value>
        struct Node;
    public:
        struct Args {
            T lr;
            T beta1;
            T beta2;
            T epsilon;
            T decay;
        };
    private:
        std::unordered_map<void*, std::any> targetNodeMap;
        Args args;
    public:
        Adam(T lr = 1E-3, T beta1 = 0.9, T beta2 = 0.999, T epsilon = 1E-8, T decay = 0);
        Adam(Args args_);
        Adam(const Adam&) = default;
        Adam(Adam&&) noexcept = default;
        ~Adam() = default;
        /* Operators */
        Adam& operator=(Adam obj) noexcept;
        /* Operations */
        void step(Diffable auto& target);
        void step(Diffable auto& target, Diffable auto&... targets);
        void clear() { targetNodeMap.clear(); }

        void swap(Adam& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] auto&& getArgs(this auto&& self) noexcept { return self.args; }
        [[nodiscard]] auto&& getLearnRate(this auto&& self) noexcept { return self.lr; }
    };

    template<Scalar T>
    Adam<T>::Adam(T lr, T beta1, T beta2, T epsilon, T decay)
            : Adam(Args{.lr = lr, .beta1 = beta1, .beta2 = beta2, .epsilon = epsilon, .decay = decay}) {}

    template<Scalar T>
    Adam<T>::Adam(Args args_) : args(args_) {
        assert(args.lr.isPositive());
        assert(args.beta1.isPositive() && args.beta1 < T(1));
        assert(args.beta2.isPositive() && args.beta2 < T(1));
        assert(args.epsilon.isPositive());
        assert(!args.decay.isNegative() && args.decay < T(1));
    }

    template<Scalar T>
    auto Adam<T>::operator=(Adam obj) noexcept -> Adam& {
        swap(obj);
        return *this;
    }

    template<Scalar T>
    void Adam<T>::step(Diffable auto& target) {
        using Target = decltype(target);
        using NodeT = Node<typename Base::template ValueT<Target>>;
        if (!args.decay.isZero())
            target.grads() += args.decay * target.values();

        void* pTarget = (void*)(&target);
        const bool exist = targetNodeMap.contains(pTarget);
        auto& any = targetNodeMap[pTarget];
        if (!exist)
            any = NodeT(args, target);

        const T beta1 = args.beta1;
        const T beta2 = args.beta2;
        auto& [m, v, beta1t, beta2t] = std::any_cast<NodeT>(any);
        m = beta1 * m + (T(1) - beta1) * target.grads();
        v = beta2 * v + (T(1) - beta2) * target.grads().squaredNorms();
        const T alpha = args.lr / (T(1) - beta1t) * sqrt(T(1) - beta2t);
        if constexpr (Vector<Target>)
            target.values() -= alpha * divide(m, sqrt(v) + args.epsilon);
        else
            target.values() -= alpha * divide(m, sqrt_elem(v) + args.epsilon);
        beta1t *= beta1;
        beta2t *= beta2;
    }

    template<Scalar T>
    void Adam<T>::step(Diffable auto& target, Diffable auto&... targets) {
        step(target);
        step(targets...);
    }

    template<Scalar T>
    void Adam<T>::swap(Adam& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        targetNodeMap.swap(obj.targetNodeMap);
        std::swap(args, obj.args);
    }

    template<Scalar T>
    template<class Value>
    struct Adam<T>::Node {
        Value m;
        Value v;
        T beta1t;
        T beta2t;

        Node() = default;
        Node(const Args& args, const auto& src);
    };

    template<Scalar T>
    template<class Value>
    Adam<T>::Node<Value>::Node(const Args& args, const auto& src) : beta1t(args.beta1), beta2t(args.beta2) {
        if constexpr (Scalar<Value>)
            m = v = 0;
        else {
            m.resize(src);
            v.resize(src);
            m.zeros();
            v.zeros();
        }
    }
}
