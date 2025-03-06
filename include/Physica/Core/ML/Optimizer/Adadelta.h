/*
 * Copyright 2025 Weibo He.
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
     * [1] arXiv:1212.5701; https://doi.org/10.48550/arXiv.1212.5701
     * [2] pytorch; https://pytorch.org/docs/stable/generated/torch.optim.Adadelta.html
     */
    template<Scalar T>
    class Adadelta {
        static_assert(!Diffable<T>);
        using This = Adadelta<T>;
    private:
        struct VectorBuffer {
            VectorND<T> v;
            VectorND<T> u;

            VectorBuffer() = default;
            VectorBuffer(size_t length) : v(length, 0), u(length, 0){}
        };

        struct MatrixBuffer {
            DenseMatrix<T, MatrixOption::Col | MatrixOption::Element> v;
            DenseMatrix<T, MatrixOption::Col | MatrixOption::Element> u;

            MatrixBuffer() = default;
            MatrixBuffer(size_t row, size_t col) : v(row, col, 0), u(row, col, 0) {}
        };

        std::unordered_map<void*, std::variant<VectorBuffer, MatrixBuffer>> targetBufferMap;
        T lr = 1;
        T rho = 0.9;
        T epsilon = 1E-6;
        T decay = 0;
    public:
        Adadelta() = default;
        Adadelta(const This&) = default;
        Adadelta(This&&) noexcept = default;
        ~Adadelta() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<Diffable U>
        void step(U& target);

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] T getLearnRate() const noexcept { return lr; }
    };

    template<Scalar T>
    template<Diffable U>
    void Adadelta<T>::step(U& target) {
        using BufferType = std::conditional<Vector<U>, VectorBuffer, MatrixBuffer>::type;
        if (!decay.isZero())
            target.grads() += decay * target.values();

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
        }

        auto& buffer = std::get<BufferType>(var);
        auto& v = buffer.v;
        auto& u = buffer.u;

        v = rho * v + (T(1) - rho) * target.grads().squaredNorms();
        if constexpr (Vector<U>)
            target.grads() = hadamard(sqrt(divide(u + epsilon, v + epsilon)), target.grads());
        else
            target.grads() = hadamard(sqrt_elem(divide(u + epsilon, v + epsilon)), target.grads());
        u = rho * u + (T(1) - rho) * target.grads().squaredNorms();
        target.values() -= lr * target.grads();
    }

    template<Scalar T>
    void Adadelta<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        targetBufferMap.swap(obj.targetBufferMap);
        lr.swap(obj.lr);
        rho.swap(obj.rho);
        epsilon.swap(obj.epsilon);
        decay.swap(obj.decay);
    }
}
