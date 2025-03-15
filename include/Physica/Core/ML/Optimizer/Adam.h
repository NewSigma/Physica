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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

namespace Physica {
    /**
     * Reference:
     * [1] arXiv:1412.6980; https://doi.org/10.48550/arXiv.1412.6980
     * [2] pytorch; https://pytorch.org/docs/stable/generated/torch.optim.Adam.html
     */
    template<Scalar T>
    class Adam {
        static_assert(!Diffable<T>);
        using This = Adam<T>;
    public:
        struct Args {
            T lr;
            T beta1;
            T beta2;
            T epsilon;
            T decay;
        };
    private:
        struct VectorBuffer {
            VectorND<T> m;
            VectorND<T> v;
            T beta1t;
            T beta2t;

            VectorBuffer() = default;
            VectorBuffer(const Args& args, size_t length) : m(length), v(length), beta1t(args.beta1), beta2t(args.beta2) {
                m.zeros();
                v.zeros();
            }
        };

        struct MatrixBuffer {
            DenseMatrix<T, MatrixOption::Col | MatrixOption::Element> m;
            DenseMatrix<T, MatrixOption::Col | MatrixOption::Element> v;
            T beta1t;
            T beta2t;

            MatrixBuffer() = default;
            MatrixBuffer(const Args& args, size_t row, size_t col) : m(row, col), v(row, col), beta1t(args.beta1), beta2t(args.beta2) {
                m.zeros();
                v.zeros();
            }
        };

        std::unordered_map<void*, std::variant<VectorBuffer, MatrixBuffer>> targetBufferMap;
        Args args;
    public:
        Adam(T lr = 1E-3, T beta1 = 0.9, T beta2 = 0.999, T epsilon = 1E-8, T decay = 0);
        Adam(Args args_);
        Adam(const Adam&) = default;
        Adam(Adam&&) noexcept = default;
        ~Adam() = default;
        /* Operators */
        Adam& operator=(Adam obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<Diffable U>
        void step(U& target);
        void clear() { targetBufferMap.clear(); }

        void swap(Adam& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const Args& getArgs() const noexcept { return args; }
        [[nodiscard]] Args& getArgs() noexcept { return args; }
        [[nodiscard]] T& getLearnRate() noexcept { return args.lr; }
        [[nodiscard]] const T& getLearnRate() const noexcept { return args.lr; }
        /* Friends */
        friend class device_obj<This>;
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
    template<Diffable U>
    void Adam<T>::step(U& target) {
        using BufferType = std::conditional<Vector<U>, VectorBuffer, MatrixBuffer>::type;
        if (!args.decay.isZero())
            target.grads() += args.decay * target.values();

        void* pTarget = (void*)(&target);
        const bool exist = targetBufferMap.count(pTarget) != 0;
        auto& var = targetBufferMap[pTarget];
        if (!exist) {
            if constexpr (Vector<U>)
                var = BufferType(args, target.getLength());
            else {
                static_assert(Matrix<U>, "[Error]: Unexpected type");
                var = BufferType(args, target.getRow(), target.getCol());
            }
        }

        const T beta1 = args.beta1;
        const T beta2 = args.beta2;
        auto& buffer = std::get<BufferType>(var);
        auto& m = buffer.m;
        auto& v = buffer.v;
        auto& beta1t = buffer.beta1t;
        auto& beta2t = buffer.beta2t;
        m = beta1 * m + (T(1) - beta1) * target.grads();
        v = beta2 * v + (T(1) - beta2) * target.grads().squaredNorms();
        const T alpha = args.lr / (T(1) - beta1t) * sqrt(T(1) - beta2t);
        if constexpr (Vector<U>)
            target.values() -= alpha * divide(m, sqrt(v) + args.epsilon);
        else
            target.values() -= alpha * divide(m, sqrt_elem(v) + args.epsilon);
        beta1t *= beta1;
        beta2t *= beta2;
    }

    template<Scalar T>
    void Adam<T>::swap(Adam& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        targetBufferMap.swap(obj.targetBufferMap);
        std::swap(args, obj.args);
    }
}
