/*
 * Copyright 2022-2025 Weibo He.
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
    template<Vector V>
    __host__ __device__ void device_obj<RValueVector<Derived>>::assign(V& target) const {
        assert_assign(target);
        if (IsHost()) {
            auto func = [source_ = asStruct(Base::getDerived()), target_ = asStruct(target)] __device__() mutable {
                const auto& source = source_.getDerived();
                auto& target = target_.getDerived();
                const size_t length = source.getLength();
                const uint32_t index = blockIdx.x * blockDim.x + threadIdx.x;
                if (index < length)
                    target[index] = source.calc(index);
            };
            CUDAExecutor::launch<MaxThreadPerBlock>(func, makeKernelConfig());
        }
        else if constexpr (IsDevice())
            assign_impl(target);
    }

    template<class Derived>
    __host__ __device__ void device_obj<RValueVector<Derived>>::assert_assign(const Vector auto& target) const noexcept {
        static_assert_assign(target);

        using V = std::remove_cvref<decltype(target)>::type;
        constexpr size_t Size1 = SizeAtCompile;
        constexpr size_t Size2 = V::SizeAtCompile;
        if constexpr (Size1 == Dynamic && Size2 == Dynamic)
            assert(getLength() == target.getLength() && "[Error]: Size mismatch between two vector");
    }

    template<class Derived>
    template<Packet Pack>
    __device__ Pack device_obj<RValueVector<Derived>>::packet(size_t index) const {
        assert(index + Pack::size() <= getLength());
        if constexpr (Scalar<Pack>)
            return calc(index);
        else
            return Pack(calc(index).value(), calc(index + 1).value());
    }

    template<class Derived>
    template<Packet Pack>
    __device__ Pack device_obj<RValueVector<Derived>>::packetPartial(size_t index, [[maybe_unused]] size_t count) const {
        assert(index + Pack::size() <= getLength());
        assert(count == 1 && "[Error]: No need to call partial version");
        if constexpr (Scalar<Pack>)
            return calc(index);
        else
            return Pack(calc(index).value(), 0_HF);
    }

    template<class Derived>
    void device_obj<RValueVector<Derived>>::reverse(const Vector auto&, const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff);
        Base::getDerived().reverse(grad);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueVector<Derived>>::transpose() const noexcept {
        return device_obj<Transpose<Derived>>(Base::getDerived());
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
            using U = std::conditional<isReverseDiff, Tv, T>::type;
            auto buffer = device_obj<VectorND<U>>(MaxThreadPerBlock);
            auto func = [v_ = asStruct(Base::getDerived()), buffer_ = asStruct(buffer)] __device__() mutable {
                const auto& v = v_.getDerived();
                auto& buffer = buffer_.getDerived();
                U local = 0;
                for (int i = threadIdx.x; i < v.getLength(); i += MaxThreadPerBlock) {
                    if constexpr (isReverseDiff)
                        local += v.calc_value(i);
                    else
                        local += v.calc(i);
                }
                buffer[threadIdx.x] = local;
            };
            CUDAExecutor::launch<MaxThreadPerBlock>(func, KernelConfig(1, MaxThreadPerBlock));

            if constexpr (IsHost()) { // To silence warnings
                CUDAExecutor::wait();
                return buffer.toHost().sum();
            }
            else
                unreachable();
        }
        else {
            auto result = T(0);
            for (size_t i = 0; i < getLength(); ++i) {
                if constexpr (isReverseDiff)
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
        if constexpr (isComplex)
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
    __device__ auto device_obj<RValueVector<Derived>>::crossProduct(const Vector auto& v) const noexcept requires(CUDA<decltype(v)>) {
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
    __device__ auto device_obj<RValueVector<Derived>>::lnSumExp(int tid, int numThread, T* __restrict shared) const -> T {
        assert(tid < numThread);
        const auto& v = Base::getDerived();
        Tv m;
        if constexpr (isComplex)
            m = values().reals().max(tid, numThread, shared);
        else
            m = values().max(tid, numThread, shared);
        return ln(exp(v - m).sum(tid, numThread, shared) + Trv(std::numeric_limits<Trv>::min())) + m; // Add min() to avoid ln(0)
    }

    template<class Derived>
    __device__ auto device_obj<RValueVector<Derived>>::mean(int tid, int numThread, T* __restrict shared) const -> T {
        return sum(tid, numThread, shared) / Trv(getLength());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueVector<Derived>>::reals() const noexcept -> RealsRtnTy {
        return RealsRtnTy(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueVector<Derived>>::imags() const noexcept {
        return device_obj<ImagVector<Derived>>(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueVector<Derived>>::squaredNorms() const noexcept {
        return device_obj<SquaredNormVector<Derived>>(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueVector<Derived>>::norms() const noexcept {
        return device_obj<NormVector<Derived>>(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueVector<Derived>>::values() const noexcept -> ValuesRtnTy {
        return ValuesRtnTy(Base::getDerived());
    }

    template<class Derived>
    template<int GradOrder>
    __host__ __device__ auto device_obj<RValueVector<Derived>>::grads() const noexcept {
        return device_obj<GradVector<Derived, GradOrder>>(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ KernelConfig device_obj<RValueVector<Derived>>::makeKernelConfig() const noexcept {
        return makeKernelConfig(getLength());
    }

    template<class Derived>
    __host__ __device__ KernelConfig device_obj<RValueVector<Derived>>::makeKernelConfig(size_t length) noexcept {
        constexpr uint32_t MaxThread = MaxThreadPerBlock;
        const uint32_t numThread = std::min<uint32_t>(length, MaxThread);
        const uint32_t numBlock = (length + numThread - 1) / numThread;
        return KernelConfig(numBlock, numThread);
    }

    template<class Derived>
    __host__ __device__ void device_obj<RValueVector<Derived>>::static_assert_assign(const Vector auto& target) noexcept {
        host_obj::static_assert_assign(target);
        static_assert(SizeAtCompile != Dynamic || CUDA<decltype(target)>, "[Error]: Host object cannot be assigned to dynamic device object");
    }

    template<class Derived>
    template<Vector V>
    __device__ void device_obj<RValueVector<Derived>>::assign_impl(V& target) const {
        if constexpr (Internal::EnableSIMD<device_obj<Derived>, V>::value) {
            constexpr static size_t Size1 = Derived::SizeAtCompile;
            constexpr static size_t Size2 = V::SizeAtCompile;
            constexpr static size_t VectorSize = Size1 > Size2 ? Size1 : Size2;
            using PacketType = device_obj<BestPacket<T, VectorSize>>::Type;
            constexpr static size_t PacketSize = PacketType::size();

            if constexpr (VectorSize != Dynamic) {
                constexpr size_t to = VectorSize / PacketSize * PacketSize;
                for (size_t i = 0; i < to; i += PacketSize)
                    target.writePacket(i, Base::getDerived().template packet<PacketType>(i));

                constexpr size_t i = VectorSize - VectorSize % PacketSize;
                if constexpr (i != VectorSize) {
                    constexpr size_t count = VectorSize - i;
                    target.writePacketPartial(i, count, Base::getDerived().template packetPartial<PacketType>(i, count));
                }
            }
            else {
                const size_t length = getLength();
                if (length == 0)
                    return;

                const size_t to = length / PacketSize * PacketSize;
                size_t i = 0;
                for (; i < to; i += PacketSize)
                    target.writePacket(i, Base::getDerived().template packet<PacketType>(i));

                if (to != length) {
                    const size_t count = length - i;
                    target.writePacketPartial(i, count, Base::getDerived().template packetPartial<PacketType>(i, count));
                }
            }
        }
        else {
            using OtherScalar = V::ScalarType;
            for (size_t i = 0; i < getLength(); ++i)
                target[i] = OtherScalar(calc(i));
        }
    }
}
