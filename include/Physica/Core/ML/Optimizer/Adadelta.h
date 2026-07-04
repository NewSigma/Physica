/*
 * Copyright 2025-2026 Weibo He.
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
     * [1] arXiv:1212.5701; https://doi.org/10.48550/arXiv.1212.5701
     * [2] pytorch; https://pytorch.org/docs/stable/generated/torch.optim.Adadelta.html
     */
    template<Scalar T>
    class Adadelta : public OptBase<T> {
        using This = Adadelta<T>;
        using Base = OptBase<T>;

        template<class Value>
        struct Node;
    private:
        std::unordered_map<void*, std::any> targetNodeMap;
        T lr = 1;
        T rho = 0.9;
        T epsilon = 1E-6;
        T decay = 0;
    public:
        Adadelta() = default;
        Adadelta(T lr_) : lr(lr_) {}
        Adadelta(const This&) = default;
        Adadelta(This&&) noexcept = default;
        ~Adadelta() = default;
        /* Operators */
        This& operator=(This obj) noexcept;
        /* Operations */
        void step(Diffable auto& target);
        void step(Diffable auto& target, Diffable auto&... targets);
        void clear() noexcept;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] auto&& getLearnRate(this auto&& self) noexcept { return self.lr; }
    };

    template<Scalar T>
    auto Adadelta<T>::operator=(Adadelta obj) noexcept -> Adadelta& {
        swap(obj);
        return *this;
    }

    template<Scalar T>
    void Adadelta<T>::step(Diffable auto& target) {
        using Target = decltype(target);
        using NodeT = Node<typename Base::template ValueT<Target>>;
        if (!decay.isZero())
            target.grads() += decay * target.values();

        void* pTarget = (void*)(&target);
        const bool exist = targetNodeMap.contains(pTarget);
        auto& any = targetNodeMap[pTarget];
        if (!exist)
            any = NodeT(target);

        auto& [v, u] = std::any_cast<NodeT>(any);
        v = rho * v + (T(1) - rho) * target.grads().squaredNorms();
        if constexpr (Vector<Target>)
            target.grads() = hadamard(sqrt(divide(u + epsilon, v + epsilon)), target.grads());
        else
            target.grads() = hadamard(sqrt_elem(divide_elem(u + epsilon, v + epsilon)), target.grads());
        u = rho * u + (T(1) - rho) * target.grads().squaredNorms();
        target.values() -= lr * target.grads();
    }

    template<Scalar T>
    void Adadelta<T>::step(Diffable auto& target, Diffable auto&... targets) {
        step(target);
        step(targets...);
    }

    template<Scalar T>
    void Adadelta<T>::clear() noexcept {
        targetNodeMap.clear();
    }

    template<Scalar T>
    void Adadelta<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        targetNodeMap.swap(obj.targetNodeMap);
        lr.swap(obj.lr);
        rho.swap(obj.rho);
        epsilon.swap(obj.epsilon);
        decay.swap(obj.decay);
    }

    template<Scalar T>
    template<class Value>
    struct Adadelta<T>::Node {
        Value v;
        Value u;

        Node() = default;
        Node(const auto& src);
    };

    template<Scalar T>
    template<class Value>
    Adadelta<T>::Node<Value>::Node(const auto& src) {
        if constexpr (Scalar<Value>)
            v = u = 0;
        else {
            v.resize(src);
            u.resize(src);
            v.zeros();
            u.zeros();
        }
    }
}
