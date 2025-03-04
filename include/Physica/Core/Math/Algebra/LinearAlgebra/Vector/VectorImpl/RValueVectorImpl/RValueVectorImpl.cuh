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
#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"
#include "Physica/PlainStruct.h"

namespace Physica {
    template<class Derived>
    template<Vector V>
    __host__ __device__ void device_obj<RValueVector<Derived>>::assign(V& target) const requires(CUDA<V>) {
        constexpr static size_t Size1 = SizeAtCompile;
        constexpr static size_t Size2 = V::SizeAtCompile;
        static_assert(Size1 == Dynamic || Size2 == Dynamic || Size1 == Size2, "[Error]: Size mismatch between two vector");
        static_assert(V::isComplex || !isComplex, "[Error]: Cannot convert a complex to a real");
        static_assert(Diffable<V> || !Diffable<This>, "[Error]: Assign a diffable vector to normal vector discards grads");
        assert(getLength() == target.getLength() && "[Error]: Size mismatch between two vector");
        if (IsHost()) {
            auto config = makeKernelConfig();
            CUDAExecutor::launch([source_ = asStruct(Base::getDerived()), target_ = asStruct(target)] __device__() mutable {
                const auto& source = source_.getDerived();
                auto& target = target_.getDerived();
                const size_t length = source.getLength();
                const uint32_t index = blockIdx.x * blockDim.x + threadIdx.x;
                if (index < length)
                    target[index] = source.calc(index);
            }, config.first, config.second);
        }
        else
            assign_impl<V>(target);
    }

    template<class Derived>
    template<Packet Pack>
    __device__ inline Pack device_obj<RValueVector<Derived>>::packet(size_t index) const {
        static_assert(Scalar<Pack>, "[Error]: Not implemented");
        return calc(index);
    }

    template<class Derived>
    template<Packet Pack>
    __device__ inline Pack device_obj<RValueVector<Derived>>::packetPartial(size_t index, size_t count) const {
        static_assert(Scalar<Pack>, "[Error]: Not implemented");
        assert(count == 1 && "[Error]: No need to call partial version");
        return calc(index);
    }

    template<class Derived>
    __host__ __device__ inline auto device_obj<RValueVector<Derived>>::transpose() const noexcept {
        return device_obj<TransposeVector<Derived>>(Base::getDerived());
    }

    template<class Derived>
    __device__ inline auto device_obj<RValueVector<Derived>>::norm() const -> Tr {
        return sqrt(Base::getDerived().squaredNorm());
    }

    template<class Derived>
    __device__ inline auto device_obj<RValueVector<Derived>>::squaredNorm() const -> Tr {
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
            const auto config = makeKernelConfig();
            device_obj<VectorND<U>> buffer(MaxThreadPerBlock);
            auto future = CUDAExecutor::launch([v_ = asStruct(Base::getDerived()), buffer_ = asStruct(buffer)] __device__() mutable {
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
            }, 1, MaxThreadPerBlock);

            if constexpr (IsHost()) { // To silence warnings
                future.wait();
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
    template<Vector V>
    __device__ inline auto device_obj<RValueVector<Derived>>::crossProduct(const device_obj<V>& v) const noexcept {
        return device_obj<CrossProduct<Derived, V>>(*this, v);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<RValueVector<Derived>>::reals() const noexcept {
        return device_obj<RealVector<Derived>>(Base::getDerived());
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
    __host__ __device__ std::pair<dim3, dim3> device_obj<RValueVector<Derived>>::makeKernelConfig() const noexcept {
        return makeKernelConfig(getLength());
    }

    template<class Derived>
    __host__ __device__ std::pair<dim3, dim3> device_obj<RValueVector<Derived>>::makeKernelConfig(size_t length) noexcept {
        constexpr uint32_t MaxThread = MaxThreadPerBlock;
        const uint32_t numThread = std::min<uint32_t>(length, MaxThread);
        const uint32_t numBlock = (length + numThread - 1) / numThread;
        return std::make_pair(dim3(numBlock), dim3(numThread));
    }

    template<class Derived>
    template<Vector V>
    __device__ void device_obj<RValueVector<Derived>>::assign_impl(V& target) const requires(CUDA<V>) {
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
