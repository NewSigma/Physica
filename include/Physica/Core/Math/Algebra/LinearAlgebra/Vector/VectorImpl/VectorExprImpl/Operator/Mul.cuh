/*
 * Copyright 2024-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/VectorExpr.cuh"
#include "Physica/Core/Parallel/ThreadBlock.cuh"

namespace Physica {
    template<Vector V, Scalar U>
    class device_obj<VectorExpr<ExprID::Mul, V, U>>
            : public device_obj<BinaryVectorExpr<ExprID::Mul, V, U>> {
        using Base = device_obj<BinaryVectorExpr<ExprID::Mul, V, U>>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using Base::isReverseDiff;
    public:
        using Base::Base;
        /* Operators */
        using Base::operator*;
        [[nodiscard]] __host__ __device__ auto operator*(Scalar auto x) const noexcept;
        [[nodiscard]] __host__ __device__ auto operator-(this auto&&) noexcept;
        /* Operations */
        __host__ __device__ void assign_add(Vector auto&& v) const;
        __device__ void assign_add(Vector auto&& v, const ThreadBlock& block) const;
        __host__ __device__ void assign_add_base(Vector auto&& v) const;
        void assign_add_cublas(Vector auto&& v) const noexcept;

        [[nodiscard]] __device__ T calc(size_t index) const;
        [[nodiscard]] __device__ Tv calc_value(size_t index) const;

        template<int Size>
        [[nodiscard]] __device__ SIMD<T, Size> packet(size_t index) const noexcept;
        template<int Size>
        [[nodiscard]] __device__ SIMD<T, Size> packet(size_t index, size_t count) const noexcept;

        [[nodiscard]] __host__ __device__ T sum() const { return getLHS().sum() * getRHS(); }
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    private:
        __host__ __device__ void assign_fma_for(Vector auto&  __restrict v) const  __restrict noexcept;
    };

    template<Vector V, Scalar U>
    __host__ __device__ auto device_obj<VectorExpr<ExprID::Mul, V, U>>::operator*(Scalar auto x) const noexcept {
        return getLHS() * (getRHS() * x);
    }

    template<Vector V, Scalar U>
    __host__ __device__ auto device_obj<VectorExpr<ExprID::Mul, V, U>>::operator-(this auto&& self) noexcept {
        using Self = decltype(self);
        return std::forward<Self>(self).getLHS() * (-std::forward<Self>(self).getRHS());
    }

    template<Vector V, Scalar U>
    __host__ __device__ void device_obj<VectorExpr<ExprID::Mul, V, U>>::assign_add(Vector auto&& v) const {
        assign_add_base(v);
    }

    template<Vector V, Scalar U>
    __device__ void device_obj<VectorExpr<ExprID::Mul, V, U>>::assign_add(Vector auto&& v, const ThreadBlock& block) const {
        size_t length = Base::getLength();
        int delta = block.getLength();
        for (int index = block.rank(); index < length; index += delta)
            v[index] += calc(index);
        block.sync();
    }

    template<Vector V, Scalar U>
    __host__ __device__ void device_obj<VectorExpr<ExprID::Mul, V, U>>::assign_add_base(Vector auto&& v) const {
        using Source = std::remove_cvref<V>::type;
        using Target = std::remove_cvref<decltype(v)>::type;
        using T1 = Source::ScalarType;
        using T2 = Target::ScalarType;
        constexpr bool LowerToFMA = std::same_as<T1, T2>;
        if constexpr (LowerToFMA) {
            v.assert_assign(Base::getDerived());
            assign_fma_for(v);
        }
        else
            Base::assign_add(v);
    }

    template<Vector V, Scalar U>
    void device_obj<VectorExpr<ExprID::Mul, V, U>>::assign_add_cublas(Vector auto&& v) const noexcept {
        using Tm = decltype(std::declval<T>().toCUDA());
        v.assert_assign(Base::getLHS());

        auto& ctx = CUDAContext::getInstance();
        ctx.setPointerMode(false);
        const size_t n = Base::getLength();
        const auto alpha = Base::getRHS().toMachine();
        const auto* x = reinterpret_cast<const Tm*>(Base::getLHS().data());
        auto* y = reinterpret_cast<Tm*>(v.data());
        if constexpr (Base::isComplex()) {
            if constexpr (T::Prec == Float32)
                check(cublasCaxpy_64(ctx, n, &alpha, x, 1, y, 1));
            else
                check(cublasZaxpy_64(ctx, n, &alpha, x, 1, y, 1));
        }
        else {
            if constexpr (T::Prec == Float32)
                check(cublasSaxpy_64(ctx, n, &alpha, x, 1, y, 1));
            else
                check(cublasDaxpy_64(ctx, n, &alpha, x, 1, y, 1));
        }
    }

    template<Vector V, Scalar U>
    __device__ auto device_obj<VectorExpr<ExprID::Mul, V, U>>::calc(size_t index) const -> T {
        if constexpr (isReverseDiff())
            return calc_value(index);
        else
            return Base::getLHS().calc(index) * Base::getRHS();
    }

    template<Vector V, Scalar U>
    __device__ auto device_obj<VectorExpr<ExprID::Mul, V, U>>::calc_value(size_t index) const -> Tv {
        return Base::getLHS().calc_value(index) * Base::getRHS().value();
    }

    template<Vector V, Scalar U>
    template<int Size>
    __device__ auto device_obj<VectorExpr<ExprID::Mul, V, U>>::packet(size_t index) const noexcept -> SIMD<T, Size> {
        return Base::getLHS().template packet<Size>(index) * SIMD<T, Size>(Base::getRHS());
    }

    template<Vector V, Scalar U>
    template<int Size>
    __device__ auto device_obj<VectorExpr<ExprID::Mul, V, U>>::packet(size_t index, size_t count) const noexcept -> SIMD<T, Size> {
           return Base::getLHS().template packet<Size>(index, count) * SIMD<T, Size>(Base::getRHS());
    }

    template<Vector V, Scalar U>
    __host__ __device__ void device_obj<VectorExpr<ExprID::Mul, V, U>>::assign_fma_for(Vector auto&  __restrict v) const  __restrict noexcept {
        if (IsHost()) {
            auto fn = [source_ = asStruct(*this), target_ = asStruct(v)] __device__() mutable {
                const auto& source = source_.getDerived();
                auto& target = target_.getDerived();
                size_t length = source.getLength();
                uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
                if (i < length) {
                    if constexpr (isReverseDiff())
                        target[i] = fma(source.getLHS().calc_value(i), Tv(source.getRHS().value()), target[i]);
                    else
                        target[i] = fma(source.getLHS().calc(i), T(source.getRHS().value()), target[i]);
                }
            };
            CUDAExecutor::launch<CUDADevAttr::DefaultThreadsPerBlock>(fn, Base::makeKernelConfig());
        }

        if constexpr (IsDevice()) {
            for (size_t i = 0; i < Base::getLength(); ++i) {
                if constexpr (isReverseDiff())
                    v[i] = fma(getLHS().calc_value(i), Tv(getRHS().value()), v[i]);
                else
                    v[i] = fma(getLHS().calc(i), T(getRHS().value()), v[i]);
            }
        }
    }
    ////////////////////////////////////////////////////////////////////
    template<Vector V1, Vector V2>
    class device_obj<VectorExpr<ExprID::Mul, V1, V2>>
            : public device_obj<BinaryVectorExpr<ExprID::Mul, V1, V2>> {
        using This = device_obj<VectorExpr<ExprID::Mul, V1, V2>>;
        using Base = device_obj<BinaryVectorExpr<ExprID::Mul, V1, V2>>;
    protected:
        using typename Base::T;
    public:
        using Base::Base;
        /* Operations */
        __host__ __device__ void assign_add(Vector auto&& v) const;

        [[nodiscard]] __device__ T calc(size_t index) const;

        template<int Size>
        [[nodiscard]] __device__ SIMD<T, Size> packet(size_t index) const noexcept;
        template<int Size>
        [[nodiscard]] __device__ SIMD<T, Size> packet(size_t index, size_t count) const noexcept;
        /* Getters */
        using Base::getLHS;
        using Base::getRHS;
    private:
        __host__ __device__ void assign_fma_for(Vector auto&  __restrict v) const  __restrict noexcept;
    };

    template<Vector V1, Vector V2>
    __host__ __device__ void device_obj<VectorExpr<ExprID::Mul, V1, V2>>::assign_add(Vector auto&& v) const {
        using Target = std::remove_cvref_t<decltype(v)>;
        using T1 = std::remove_cvref_t<V1>::ScalarType;
        using T2 = std::remove_cvref_t<V2>::ScalarType;
        using U = Target::ScalarType;
        constexpr bool LowerToFMA = std::same_as<T1, T2> && std::same_as<T, U>;
        if constexpr (LowerToFMA) {
            v.assert_assign(Base::getDerived());
            assign_fma_for(v);
        }
        else
            Base::assign_add(v);
    }

    template<Vector V1, Vector V2>
    __device__ auto device_obj<VectorExpr<ExprID::Mul, V1, V2>>::calc(size_t index) const -> T {
        return Base::getLHS().calc(index) * Base::getRHS().calc(index);
    }

    template<Vector V1, Vector V2>
    template<int Size>
    __device__ auto device_obj<VectorExpr<ExprID::Mul, V1, V2>>::packet(size_t index) const noexcept -> SIMD<T, Size> {
        return Base::getLHS().template packet<Size>(index)
             * Base::getRHS().template packet<Size>(index);
    }

    template<Vector V1, Vector V2>
    template<int Size>
    __device__ auto device_obj<VectorExpr<ExprID::Mul, V1, V2>>::packet(size_t index, size_t count) const noexcept -> SIMD<T, Size> {
        return Base::getLHS().template packet<Size>(index, count)
             * Base::getRHS().template packet<Size>(index, count);
    }

    template<Vector V1, Vector V2>
    __host__ __device__ void device_obj<VectorExpr<ExprID::Mul, V1, V2>>::assign_fma_for(Vector auto&  __restrict v) const  __restrict noexcept {
        if (IsHost()) {
            auto fn = [source_ = asStruct(*this), target_ = asStruct(v)] __device__() mutable {
                const auto& source = source_.getDerived();
                auto& target = target_.getDerived();
                size_t length = source.getLength();
                uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
                if (i < length) {
                    if constexpr (Base::isReverseDiff())
                        target[i] = fma(source.getLHS().calc_value(i), Tv(source.getRHS().calc_value(i)), target[i]);
                    else
                        target[i] = fma(source.getLHS().calc(i), T(source.getRHS().calc(i)), target[i]);
                }
            };
            CUDAExecutor::launch<CUDADevAttr::DefaultThreadsPerBlock>(fn, Base::makeKernelConfig());
        }

        if constexpr (IsDevice()) {
            for (size_t i = 0; i < Base::getLength(); ++i) {
                if constexpr (Base::isReverseDiff())
                    v[i] = fma(getLHS().calc_value(i), Tv(getRHS().calc_value(i)), v[i]);
                else
                    v[i] = fma(getLHS().calc(i), T(getRHS().calc(i)), v[i]);
            }
        }
    }
    ////////////////////////////////////////////////////////////////////
    template<Vector V, Scalar U>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator*(V&& v, U&& x) noexcept requires(DeviceObj<V>) {
        return device_obj<VectorExpr<ExprID::Mul, remove_device_obj_t<V&&>, U&&>>(std::forward<V>(v), std::forward<U>(x));
    }

    template<Scalar U, Vector V>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto operator*(U&& x, V&& v) noexcept requires(DeviceObj<V>) {
        return std::forward<V>(v) * std::forward<U>(x);
    }

    template<Vector V1, Vector V2>
    [[nodiscard, gnu::always_inline]] __host__ __device__ auto hadamard(V1&& v1, V2&& v2) noexcept requires(DeviceObj<V1> && DeviceObj<V2>) {
        if constexpr (!canonicalized(v1, v2))
            return hadamard(std::forward<V2>(v2), std::forward<V1>(v1));
        else
            return device_obj<VectorExpr<ExprID::Mul, remove_device_obj_t<V1&&>, remove_device_obj_t<V2&&>>>(std::forward<V1>(v1), std::forward<V2>(v2));
    }
}
