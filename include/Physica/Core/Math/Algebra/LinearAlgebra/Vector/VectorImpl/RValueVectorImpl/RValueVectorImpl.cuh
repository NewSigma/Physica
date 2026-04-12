/*
 * Copyright 2022-2026 Weibo He.
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

#include "../RValueVector.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/PlainStruct.h"

namespace Physica {
    template<class Derived>
    __host__ __device__ auto device_obj<RValueVector<Derived>>::operator*(this auto&& self, Matrix auto&& m) noexcept {
        using Self = decltype(self);
        using M = decltype(m);
        static_assert(is_device_obj<M>::value, "[Error]: host-device mismatch");
        return device_obj<GEVM<Self&&, M&&>>(std::forward<Self>(self), std::forward<M>(m));
    }

    template<class Derived>
    __host__ __device__ void device_obj<RValueVector<Derived>>::assign(Vector auto&& target) const {
        target.assert_assign(Base::getDerived());
        if (IsHost()) {
            auto func = [source_ = asStruct(Base::getDerived()), target_ = asStruct(target)] __device__() mutable {
                const auto& source = source_.getDerived();
                auto& target = target_.getDerived();
                size_t length = source.getLength();
                uint32_t index = blockIdx.x * blockDim.x + threadIdx.x;
                if (index < length)
                    target[index] = source.calc(index);
            };
            CUDAExecutor::launch<CUDADevAttr::DefaultThreadsPerBlock>(func, makeKernelConfig());
        }
        else if constexpr (IsDevice())
            assign_impl(target);
    }

    template<class Derived>
    __device__ void device_obj<RValueVector<Derived>>::assign(Vector auto&& target, const ThreadBlock& block) const {
        for (int i = block.rank(); i < getLength(); i += block.getLength())
            target[i] = calc(i);
        block.sync();
    }

    template<class Derived>
    __host__ __device__ void device_obj<RValueVector<Derived>>::assign_add(Vector auto& target) const {
        target.assert_assign(Base::getDerived());
        if (IsHost()) {
            auto func = [source_ = asStruct(Base::getDerived()), target_ = asStruct(target)] __device__() mutable {
                const auto& source = source_.getDerived();
                auto& target = target_.getDerived();
                size_t length = source.getLength();
                uint32_t index = blockIdx.x * blockDim.x + threadIdx.x;
                if (index < length)
                    target[index] += source.calc(index);
            };
            CUDAExecutor::launch<CUDADevAttr::DefaultThreadsPerBlock>(func, makeKernelConfig());
        }
        else if constexpr (IsDevice())
            noImpl();
    }

    template<class Derived>
    __host__ __device__ void device_obj<RValueVector<Derived>>::assign_add_base(Vector auto& target) const {
        assign_add(target);
    }

    template<class Derived>
    __host__ __device__ void device_obj<RValueVector<Derived>>::assert_assign(const Vector auto& source) const noexcept {
        static_assert_assign(source);

        using V = std::remove_cvref<decltype(source)>::type;
        constexpr size_t Size1 = SizeAtCompile;
        constexpr size_t Size2 = V::SizeAtCompile;
        if constexpr (Size1 == Dynamic && Size2 == Dynamic)
            assert(getLength() == source.getLength() && "[Error]: Size mismatch between two vector");
    }

    template<class Derived>
    __device__ auto device_obj<RValueVector<Derived>>::calc(size_t index) const -> T {
        return Base::getDerived().calc(index);
    }

    template<class Derived>
    __device__ auto device_obj<RValueVector<Derived>>::calc_value(size_t index) const -> Tv {
        return Base::getDerived().values().calc(index);
    }

    template<class Derived>
    template<int Size>
    __device__ auto device_obj<RValueVector<Derived>>::packet(size_t index) const noexcept -> SIMD<T, Size> {
        assert(index + Size <= getLength());
        return SIMD<T, Size>(calc(index).value(), calc(index + 1).value());
    }

    template<class Derived>
    template<int Size>
    __device__ auto device_obj<RValueVector<Derived>>::packet(size_t index, [[maybe_unused]] size_t count) const noexcept -> SIMD<T, Size> {
        assert(index + Size <= getLength());
        assert(count == 1 && "[Error]: No need to call partial version");
        return SIMD<T, Size>(calc(index).value(), 0_HF);
    }

    template<class Derived>
    void device_obj<RValueVector<Derived>>::reverse(const Vector auto&, const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff());
        Base::getDerived().reverse(grad);
    }

    template<class Derived>
    __host__ __device__ decltype(auto) device_obj<RValueVector<Derived>>::conjugate(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (isComplex()) {
            using V = remove_device_obj<Self>::type;
            return device_obj<Conjugate<V>>(std::forward<Self>(self));
        }
        else
            return std::forward<Self>(self);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueVector<Derived>>::transpose(this auto&& self) noexcept {
        using Self = decltype(self);
        using V = remove_device_obj<Self>::type;
        return device_obj<Transpose<V>>(std::forward<Self>(self));
    }

    template<class Derived>
    __device__ auto device_obj<RValueVector<Derived>>::norm() const -> Tr {
        return sqrt(Base::getDerived().squaredNorm());
    }

    template<class Derived>
    __device__ auto device_obj<RValueVector<Derived>>::squaredNorm() const -> Tr {
        auto result = Tr(0);
        for (size_t i = 0; i < getLength(); ++i)
            result += calc(i).squaredNorm();
        return result;
    }

    template<class Derived>
    __device__ auto device_obj<RValueVector<Derived>>::max() const -> T {
        assert(getLength() != 0);
        T result = calc(0);
        for (size_t i = 1; i < getLength(); ++i) {
            T temp = calc(i);
            if (result < temp)
                result = temp;
        }
        return result;
    }

    template<class Derived>
    __device__ auto device_obj<RValueVector<Derived>>::min() const -> T {
        assert(getLength() != 0);
        T result = calc(0);
        for (size_t i = 1; i < getLength(); ++i) {
            T temp = calc(i);
            if (result > temp)
                result = temp;
        }
        return result;
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueVector<Derived>>::sum() const -> T {
        assert(getLength() != 0);
        if (IsHost()) {
            using U = std::conditional<isReverseDiff(), Tv, T>::type;
            auto numThreads = std::min<size_t>(getLength(), CUDADevAttr::DefaultThreadsPerBlock);
            auto buffer = device_obj<VectorND<U>>(numThreads);
            auto func = [v_ = asStruct(Base::getDerived()), buffer_ = asStruct(buffer)] __device__() mutable {
                const auto& v = v_.getDerived();
                auto& buffer = buffer_.getDerived();
                U local = 0;
                for (int i = (int)threadIdx.x; i < v.getLength(); i += (int)blockDim.x) {
                    if constexpr (isReverseDiff())
                        local += v.calc_value(i);
                    else
                        local += v.calc(i);
                }
                buffer[threadIdx.x] = local;
            };
            CUDAExecutor::launch<CUDADevAttr::DefaultThreadsPerBlock>(func, KernelConfig(1, numThreads));

            if constexpr (IsHost()) { // To silence warnings
                CUDAExecutor::wait();
                return buffer.toHost().sum();
            }
            else
                unreachable();
        }
        else if constexpr (IsDevice()) {
            auto result = T(0);
            for (size_t i = 0; i < getLength(); ++i) {
                if constexpr (isReverseDiff())
                    result += calc_value(i);
                else
                    result += calc(i);
            }
            return result;
        }
    }

    template<class Derived>
    __device__ auto device_obj<RValueVector<Derived>>::mean() const -> T {
        return sum() / Trv(getLength());
    }

    template<class Derived>
    __device__ auto device_obj<RValueVector<Derived>>::variance() const -> T {
        const auto& x = Base::getDerived();
        const size_t length = getLength();
        const auto x_mean = mean();
        const auto expr = x - x_mean;
        const auto expr2 = square(expr);
        return expr2.sum() / Trv(length - 1);
    }

    template<class Derived>
    __device__ auto device_obj<RValueVector<Derived>>::variance(const T& prior_mean) const -> T {
        const auto& x = Base::getDerived();
        const size_t length = getLength();
        return (x - prior_mean).squaredNorm() / Trv(length - 1);
    }

    template<class Derived>
    __device__ auto device_obj<RValueVector<Derived>>::deviation() const -> T {
        return sqrt(variance());
    }

    template<class Derived>
    __device__ auto device_obj<RValueVector<Derived>>::deviation(const T& prior_mean) const -> T {
        return sqrt(variance(prior_mean));
    }

    template<class Derived>
    __device__ auto device_obj<RValueVector<Derived>>::lnSumExp() const -> T {
        const auto& v = Base::getDerived();
        Tv m;
        if constexpr (isComplex())
            m = values().reals().max();
        else
            m = values().max();
        return ln(exp(v - m).sum() + Trv(std::numeric_limits<Trv>::min())) + m; // Add min() to avoid ln(0)
    }

    template<class Derived>
    __device__ auto device_obj<RValueVector<Derived>>::crossEntropy(size_t index) const -> T {
        assert(index < getLength() && "[Error]: Index overflow");
        const auto& v = Base::getDerived();
        return (v - calc(index)).lnSumExp();
    }

    template<class Derived>
    __device__ auto device_obj<RValueVector<Derived>>::lnSoftmax(size_t index) const -> T {
        assert(index < getLength() && "[Error]: Index overflow");
        return -crossEntropy(index);
    }

    template<class Derived>
    __device__ auto device_obj<RValueVector<Derived>>::softmax(size_t index) const -> T {
        assert(index < getLength() && "[Error]: Index overflow");
        return exp(lnSoftmax(index));
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueVector<Derived>>::prod() const noexcept -> T {
        assert(getLength() != 0);
        if (IsHost()) {
            using U = std::conditional<isReverseDiff(), Tv, T>::type;
            auto numThreads = std::min<size_t>(getLength(), CUDADevAttr::DefaultThreadsPerBlock);
            auto buffer = device_obj<VectorND<U>>(numThreads);
            auto func = [v_ = asStruct(Base::getDerived()), buffer_ = asStruct(buffer)] __device__() mutable {
                const auto& v = v_.getDerived();
                auto& buffer = buffer_.getDerived();
                U local = 1;
                for (int i = (int)threadIdx.x; i < v.getLength(); i += (int)blockDim.x) {
                    if constexpr (isReverseDiff())
                        local *= v.calc_value(i);
                    else
                        local *= v.calc(i);
                }
                buffer[threadIdx.x] = local;
            };
            CUDAExecutor::launch<CUDADevAttr::DefaultThreadsPerBlock>(func, KernelConfig(1, numThreads));

            if constexpr (IsHost()) { // To silence warnings
                CUDAExecutor::wait();
                return buffer.toHost().prod();
            }
            else
                unreachable();
        }
        else if constexpr (IsDevice()) {
            T result = calc(0);
            for(size_t i = 1; i < getLength(); ++i)
                result *= calc(i);
            return result;
        }
    }

    template<class Derived>
    __device__ auto device_obj<RValueVector<Derived>>::cross(const Vector auto& v) const noexcept {
        static_assert(is_device_obj_v<decltype(v)>, "[Error]: host-device mismatch");
        return device_obj<CrossProduct<Derived, decltype(v)>>(*this, v);
    }

    template<class Derived>
    __device__ auto device_obj<RValueVector<Derived>>::max(int tid, int numThread, T* __restrict shared) const -> T {
        T local = std::numeric_limits<T>::lowest();
        for (int i = tid; i < getLength(); i += numThread)
            local = std::max(local, calc(i));
        shared[tid] = local;
        __syncthreads();

        for (int i = (numThread + 1) / 2; i > 0; i /= 2) {
            if ((tid < i) && (tid + i < numThread))
                shared[tid] = std::max(shared[tid], shared[tid + i]);
            __syncthreads();
        }
        return shared[0];
    }

    template<class Derived>
    __device__ auto device_obj<RValueVector<Derived>>::sum(int tid, int numThread, T* __restrict shared) const -> T {
        T local = 0;
        for (int i = tid; i < getLength(); i += numThread)
            local += calc(i);
        shared[tid] = local;
        __syncthreads();

        for (int i = (numThread + 1) / 2; i > 0; i /= 2) {
            if ((tid < i) && (tid + i < numThread))
                shared[tid] += shared[tid + i];
            __syncthreads();
        }
        return shared[0];
    }

    template<class Derived>
    __device__ auto device_obj<RValueVector<Derived>>::mean(int tid, int numThread, T* __restrict shared) const -> T {
        return sum(tid, numThread, shared) / Trv(getLength());
    }

    template<class Derived>
    __device__ auto device_obj<RValueVector<Derived>>::lnSumExp(int tid, int numThread, T* __restrict shared) const -> T {
        assert(tid < numThread);
        const auto& v = Base::getDerived();
        Tv m;
        if constexpr (isComplex())
            m = values().reals().max(tid, numThread, shared);
        else
            m = values().max(tid, numThread, shared);
        return ln(exp(v - m).sum(tid, numThread, shared) + Trv(std::numeric_limits<Trv>::min())) + m; // Add min() to avoid ln(0)
    }

    template<class Derived>
    __host__ __device__ decltype(auto) device_obj<RValueVector<Derived>>::reals(this auto&& self) noexcept {
        using Self = decltype(self);
        using V = remove_device_obj<Self>::type;
        if constexpr (isComplex())
            return device_obj<RealVectorR<V>>(std::forward<Self>(self));
        else
            return std::forward<Self>(self);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueVector<Derived>>::imags(this auto&& self) noexcept {
        using Self = decltype(self);
        using V = remove_device_obj<Self>::type;
        return device_obj<ImagVector<V>>(std::forward<Self>(self));
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueVector<Derived>>::squaredNorms(this auto&& self) noexcept {
        using Self = decltype(self);
        using V = remove_device_obj<Self>::type;
        return device_obj<SquaredNormVector<V>>(std::forward<Self>(self));
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueVector<Derived>>::norms(this auto&& self) noexcept {
        using Self = decltype(self);
        using V = remove_device_obj<Self>::type;
        return device_obj<NormVector<V>>(std::forward<Self>(self));
    }

    template<class Derived>
    __host__ __device__ decltype(auto) device_obj<RValueVector<Derived>>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        using V = remove_device_obj<Self>::type;
        if constexpr (isDiffable())
            return device_obj<ValueVector<V>>(std::forward<Self>(self));
        else
            return std::forward<Self>(self);
    }

    template<class Derived>
    template<int GradOrder>
    __host__ __device__ decltype(auto) device_obj<RValueVector<Derived>>::grads(this auto&& self) noexcept {
        using Self = decltype(self);
        using V = remove_device_obj<Self>::type;
        return device_obj<GradVector<V, GradOrder>>(std::forward<Self>(self));
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<RValueVector<Derived>>::isComplex() noexcept {
        return ScalarType::isComplex();
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<RValueVector<Derived>>::isDiffable() noexcept {
        return ScalarType::isDiffable();
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<RValueVector<Derived>>::isForwardDiff() noexcept {
        return ScalarType::isForwardDiff();
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<RValueVector<Derived>>::isReverseDiff() noexcept {
        return ScalarType::isReverseDiff();
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<RValueVector<Derived>>::isLValueVector() noexcept {
        return false;
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<RValueVector<Derived>>::isCompact() noexcept {
        return host_obj::isCompact();
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<RValueVector<Derived>>::isSparse() noexcept {
        return host_obj::isSparse();
    }

    template<class Derived>
    __host__ __device__ consteval bool device_obj<RValueVector<Derived>>::isFastPacket() noexcept {
        return host_obj::isFastPacket();
    }

    template<class Derived>
    __host__ __device__ consteval size_t device_obj<RValueVector<Derived>>::getSizeAtCompile() noexcept {
        return host_obj::getSizeAtCompile();
    }

    template<class Derived>
    __host__ __device__ consteval size_t device_obj<RValueVector<Derived>>::getSizeAtCompile(const Vector auto& hint) noexcept {
        return host_obj::getSizeAtCompile(hint);
    }

    template<class Derived>
    __host__ __device__ KernelConfig device_obj<RValueVector<Derived>>::makeKernelConfig() const noexcept {
        return makeKernelConfig(getLength());
    }

    template<class Derived>
    __host__ __device__ KernelConfig device_obj<RValueVector<Derived>>::makeKernelConfig(size_t length) noexcept {
        assert(length > 0 && "[Error]: Do not schedule empty work");
        const uint32_t numThread = std::min<uint32_t>(length, CUDADevAttr::DefaultThreadsPerBlock);
        const uint32_t numBlock = (length + numThread - 1) / numThread;
        return KernelConfig(numBlock, numThread);
    }
    // Redeclare to expose it to base classes
    template<class Derived>
    __host__ __device__ consteval void device_obj<RValueVector<Derived>>::static_assert_assign(const Scalar auto& source) noexcept {
        host_obj::static_assert_assign(source);
    }

    template<class Derived>
    __host__ __device__ consteval void device_obj<RValueVector<Derived>>::static_assert_assign(const Vector auto& source) noexcept {
        static_assert(SizeAtCompile != Dynamic || DeviceObj<decltype(source)>, "[Error]: Host object cannot be assigned to dynamic device object");
        host_obj::static_assert_assign(source);
    }

    template<class Derived>
    template<Vector V>
    __device__ void device_obj<RValueVector<Derived>>::assign_impl(V& target) const {
        const auto& v0 = Base::getDerived();
        if constexpr (Internal::EnableSIMD<device_obj<Derived>, V>::value) {
            constexpr size_t Length = std::max(Derived::SizeAtCompile, V::SizeAtCompile);
            constexpr int Size = device_obj<BestPacket<T, Length>>::Size;
            if constexpr (Length != Dynamic) {
                constexpr size_t to = Length / Size * Size;
                for (size_t i = 0; i < to; i += Size)
                    target.writePacket(v0.template packet<Size>(i), i);

                for (size_t i = Length - Length % Size; i < Length; ++i)
                    target[i] = v0.calc(i);
            }
            else {
                const size_t length = getLength();
                const size_t to = length / Size * Size;
                size_t i = 0;
                for (; i < to; i += Size)
                    target.writePacket(v0.template packet<Size>(i), i);

                for (; i < length; ++i)
                    target[i] = v0.calc(i);
            }
        }
        else {
            using OtherScalar = V::ScalarType;
            for (size_t i = 0; i < getLength(); ++i)
                target[i] = OtherScalar(calc(i));
        }
    }
}
