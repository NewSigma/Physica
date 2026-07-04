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
    template<class Derived, Scalar ScalarT>
    __host__ __device__ auto device_obj<RValueVector<Derived, ScalarT>>::operator*(this auto&& self, Scalar auto&& x) noexcept {
        using V = decltype(self);
        using U = decltype(x);
        return device_obj<VectorExpr<ExprID::Mul, remove_device_obj_t<V>, U>>(std::forward<V>(self), std::forward<U>(x));
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ auto device_obj<RValueVector<Derived, ScalarT>>::operator*(this auto&& self, Matrix auto&& m) noexcept {
        using Self = decltype(self);
        using M = decltype(m);
        static_assert(is_device_obj<M>::value, "[Error]: host-device mismatch");
        return device_obj<GEVM<Self&&, M&&>>(std::forward<Self>(self), std::forward<M>(m));
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ auto device_obj<RValueVector<Derived, ScalarT>>::operator-(this auto&& self) noexcept {
        using Self = decltype(self);
        return device_obj<VectorExpr<ExprID::Minus, remove_device_obj_t<Self>>>(std::forward<Self>(self));
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ void device_obj<RValueVector<Derived, ScalarT>>::assign(this const auto& self, Vector auto&& target) {
        target.assert_assign(self);
        if (IsHost()) {
            auto func = [source_ = asStruct(self), target_ = asStruct(target)] __device__() mutable {
                const auto& source = source_.getDerived();
                auto& target = target_.getDerived();
                size_t i = blockIdx.x * blockDim.x + threadIdx.x;
                if (i < source.getLength())
                    target[i] = source.calc(i);
            };
            CUDAExecutor::launch<CUDADevAttr::DefaultThreadsPerBlock>(func, self.makeKernelConfig());
        }

        if constexpr (IsDevice())
            self.assign(target, ThreadBlock<1>{});
    }

    template<class Derived, Scalar ScalarT>
    __device__ void device_obj<RValueVector<Derived, ScalarT>>::assign(this const auto& self, Vector auto&& target, instanceof_x<ThreadBlock> auto block) {
        target.assert_assign(self);
        if constexpr (Internal::EnableSIMD<device_obj<Derived>, decltype(target)>::value) {
            constexpr size_t Length = getSizeAtCompile(target);
            constexpr int Size = device_obj<BestPacket<T, Length>>::Size;
            if constexpr (Length != Dynamic) {
                constexpr size_t to = Length / Size * Size;
                for (size_t i = 0; i < to; i += Size)
                    target.writePacket(self.template packet<Size>(i), i);

                for (size_t i = Length - Length % Size; i < Length; ++i)
                    target[i] = self.calc(i);
            }
            else {
                const size_t length = self.getLength();
                const size_t to = length / Size * Size;
                size_t i = 0;
                for (; i < to; i += Size)
                    target.writePacket(self.template packet<Size>(i), i);

                for (; i < length; ++i)
                    target[i] = self.calc(i);
            }
        }
        else {
            const size_t length = self.getLength();
            for (size_t i = block.tid(); i < length; i += block.getNumThread())
                target[i] = self.calc(i);
            block.sync();
        }
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ void device_obj<RValueVector<Derived, ScalarT>>::assign_add(Vector auto& target) const {
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

        if constexpr (IsDevice())
            Base::getDerived().assign_add(target, ThreadBlock<1>{});
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ void device_obj<RValueVector<Derived, ScalarT>>::assign_add_base(Vector auto& target) const {
        assign_add(target);
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ void device_obj<RValueVector<Derived, ScalarT>>::assert_assign(const Vector auto& source) const noexcept {
        static_assert_assign(source);
        if constexpr (std::same_as<device_obj<Derived>, std::remove_cvref_t<decltype(source)>>)
            assert(this != &source && "[Error]: Self assign is likely a bug");
        if constexpr (getSizeAtCompile() == Dynamic && source.getSizeAtCompile() == Dynamic)
            assert(getLength() == source.getLength() && "[Error]: Size mismatch between two vector");
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ constexpr KernelConfig device_obj<RValueVector<Derived, ScalarT>>::makeKernelConfig() const noexcept {
        const size_t length = getLength();
        const uint32_t numThread = std::min<uint32_t>(length, CUDADevAttr::DefaultThreadsPerBlock);
        const uint32_t numBlock = (length + numThread - 1) / numThread;
        return KernelConfig (numBlock, numThread);
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::calc(size_t index) const -> T {
        return calc(index, ThreadBlock<1>{});
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::calc(size_t index, instanceof_x<ThreadBlock> auto block) const -> T {
        return Base::getDerived().calc(index, block);
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::calc_value(size_t index) const -> Tv {
        return calc_value(index, ThreadBlock<1>{});
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::calc_value(size_t index, instanceof_x<ThreadBlock> auto block) const -> Tv {
        return Base::getDerived().values().calc(index, block);
    }

    template<class Derived, Scalar ScalarT>
    template<int Size>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::packet(this const auto& self, size_t index) noexcept -> SIMD<T, Size> {
        assert(index + Size <= self.getLength());
        return SIMD<T, Size>(self.calc(index).value(), self.calc(index + 1).value());
    }

    template<class Derived, Scalar ScalarT>
    template<int Size>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::packet(this const auto& self, size_t index, [[maybe_unused]] size_t count) noexcept -> SIMD<T, Size> {
        assert(index + Size <= self.getLength());
        assert(count == 1 && "[Error]: No need to call partial version");
        return SIMD<T, Size>(self.calc(index).value(), 0_HF);
    }

    template<class Derived, Scalar ScalarT>
    void device_obj<RValueVector<Derived, ScalarT>>::reverse(const Vector auto&, const Vector auto& grad) const noexcept {
        static_assert(isReverseDiff());
        Base::getDerived().reverse(grad);
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ decltype(auto) device_obj<RValueVector<Derived, ScalarT>>::conjugate(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (isComplex()) {
            using V = remove_device_obj<Self>::type;
            return device_obj<Conjugate<V>>(std::forward<Self>(self));
        }
        else
            return std::forward<Self>(self);
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ auto device_obj<RValueVector<Derived, ScalarT>>::transpose(this auto&& self) noexcept {
        using Self = decltype(self);
        using V = remove_device_obj<Self>::type;
        return device_obj<Transpose<V>>(std::forward<Self>(self));
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ auto device_obj<RValueVector<Derived, ScalarT>>::hermite(this auto&& self) noexcept {
        using Self = decltype(self);
        if constexpr (isComplex()) {
            using V = remove_device_obj<Self>::type;
            return device_obj<Hermite<V>>(std::forward<Self>(self));
        }
        else
            return std::forward<Self>(self).transpose();
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::norm() const -> Tr {
        return sqrt(Base::getDerived().squaredNorm());
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::squaredNorm(this const auto& self) -> Tr {
        auto result = Tr(0);
        for (size_t i = 0; i < self.getLength(); ++i)
            result += self.calc(i).squaredNorm();
        return result;
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::max() const -> T {
        return max(ThreadBlock<1>{});
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::min() const -> T {
        return min(ThreadBlock<1>{});
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ auto device_obj<RValueVector<Derived, ScalarT>>::sum() const -> T {
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

        if constexpr (IsDevice())
            return sum(ThreadBlock<1>{});
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::mean() const -> T {
        return sum() / Trv(getLength());
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::variance() const -> T {
        const auto& x = Base::getDerived();
        const size_t length = getLength();
        const auto x_mean = mean();
        const auto expr = x - x_mean;
        const auto expr2 = square(expr);
        return expr2.sum() / Trv(length - 1);
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::variance(const T& prior_mean) const -> T {
        const auto& x = Base::getDerived();
        const size_t length = getLength();
        return (x - prior_mean).squaredNorm() / Trv(length - 1);
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::deviation() const -> T {
        return sqrt(variance());
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::deviation(const T& prior_mean) const -> T {
        return sqrt(variance(prior_mean));
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::lnSumExp() const -> T {
        const auto& v = Base::getDerived();
        Tv m;
        if constexpr (isComplex())
            m = values().reals().max();
        else
            m = values().max();
        return ln(exp(v - m).sum() + Trv(std::numeric_limits<Trv>::min())) + m; // Add min() to avoid ln(0)
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::crossEntropy(this const auto& self, size_t index) -> T {
        assert(index < self.getLength() && "[Error]: Index overflow");
        return (self - self.calc(index)).lnSumExp();
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::lnSoftmax(size_t index) const -> T {
        assert(index < getLength() && "[Error]: Index overflow");
        return -crossEntropy(index);
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::softmax(size_t index) const -> T {
        assert(index < getLength() && "[Error]: Index overflow");
        return exp(lnSoftmax(index));
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ auto device_obj<RValueVector<Derived, ScalarT>>::prod(this const auto& self) noexcept -> T {
        assert(self.getLength() != 0);
        if (IsHost()) {
            using U = std::conditional<isReverseDiff(), Tv, T>::type;
            auto numThreads = std::min<size_t>(self.getLength(), CUDADevAttr::DefaultThreadsPerBlock);
            auto buffer = device_obj<VectorND<U>>(numThreads);
            auto func = [v_ = asStruct(self), buffer_ = asStruct(buffer)] __device__() mutable {
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

        if constexpr (IsDevice()) {
            T result = self.calc(0);
            for(size_t i = 1; i < self.getLength(); ++i)
                result *= self.calc(i);
            return result;
        }
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::max(this const auto& self, instanceof_x<ThreadBlock> auto block) -> T {
        assert(self.getLength() != 0);
        const int numThread = block.getNumThread();
        T local = std::numeric_limits<T>::lowest();
        for (int i = block.tid(); i < self.getLength(); i += numThread)
            local = std::max(local, self.calc(i));
        return block.max(local);
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::min(this const auto& self, instanceof_x<ThreadBlock> auto block) -> T {
        assert(self.getLength() != 0);
        const int numThread = block.getNumThread();
        T local = std::numeric_limits<T>::max();
        for (int i = block.tid(); i < self.getLength(); i += numThread)
            local = std::min(local, self.calc(i));
        return block.min(local);
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::sum(this const auto& self, instanceof_x<ThreadBlock> auto block) -> T {
        assert(self.getLength() != 0);
        const int numThread = block.getNumThread();
        T local = 0;
        for (int i = block.tid(); i < self.getLength(); i += numThread)
            local += self.calc(i);
        return block.sum(local);
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::mean(instanceof_x<ThreadBlock> auto block) const -> T {
        return sum(block) / Trv(getLength());
    }

    template<class Derived, Scalar ScalarT>
    __device__ auto device_obj<RValueVector<Derived, ScalarT>>::lnSumExp(instanceof_x<ThreadBlock> auto block) const -> T {
        Tv m;
        if constexpr (isComplex())
            m = values().reals().max(block);
        else
            m = values().max(block);
        return ln(exp(Base::getDerived() - m).sum(block) + Trv(std::numeric_limits<Trv>::min())) + m; // Add min() to avoid ln(0)
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ decltype(auto) device_obj<RValueVector<Derived, ScalarT>>::reals(this auto&& self) noexcept {
        using Self = decltype(self);
        using V = remove_device_obj<Self>::type;
        if constexpr (isComplex())
            return device_obj<RealVectorR<V>>(std::forward<Self>(self));
        else
            return std::forward<Self>(self);
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ auto device_obj<RValueVector<Derived, ScalarT>>::imags(this auto&& self) noexcept {
        using Self = decltype(self);
        using V = remove_device_obj<Self>::type;
        return device_obj<ImagVectorR<V>>(std::forward<Self>(self));
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ auto device_obj<RValueVector<Derived, ScalarT>>::squaredNorms(this auto&& self) noexcept {
        using Self = decltype(self);
        using V = remove_device_obj<Self>::type;
        return device_obj<SquaredNormVector<V>>(std::forward<Self>(self));
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ auto device_obj<RValueVector<Derived, ScalarT>>::norms(this auto&& self) noexcept {
        using Self = decltype(self);
        using V = remove_device_obj<Self>::type;
        return device_obj<NormVector<V>>(std::forward<Self>(self));
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ decltype(auto) device_obj<RValueVector<Derived, ScalarT>>::values(this auto&& self) noexcept {
        using Self = decltype(self);
        using V = remove_device_obj<Self>::type;
        if constexpr (isDiffable())
            return device_obj<ValueVector<V>>(std::forward<Self>(self));
        else
            return std::forward<Self>(self);
    }

    template<class Derived, Scalar ScalarT>
    template<int GradOrder>
    __host__ __device__ decltype(auto) device_obj<RValueVector<Derived, ScalarT>>::grads(this auto&& self) noexcept {
        using Self = decltype(self);
        using V = remove_device_obj<Self>::type;
        return device_obj<GradVector<V, GradOrder>>(std::forward<Self>(self));
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool device_obj<RValueVector<Derived, ScalarT>>::isComplex() noexcept {
        return ScalarType::isComplex();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool device_obj<RValueVector<Derived, ScalarT>>::isDiffable() noexcept {
        return ScalarType::isDiffable();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool device_obj<RValueVector<Derived, ScalarT>>::isForwardDiff() noexcept {
        return ScalarType::isForwardDiff();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool device_obj<RValueVector<Derived, ScalarT>>::isReverseDiff() noexcept {
        return ScalarType::isReverseDiff();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool device_obj<RValueVector<Derived, ScalarT>>::isLValueVector() noexcept {
        return false;
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool device_obj<RValueVector<Derived, ScalarT>>::isCompact() noexcept {
        return Derived::isCompact();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool device_obj<RValueVector<Derived, ScalarT>>::isSparse() noexcept {
        return Derived::isSparse();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool device_obj<RValueVector<Derived, ScalarT>>::isFastAssign() noexcept {
        return Derived::isFastAssign();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool device_obj<RValueVector<Derived, ScalarT>>::isFastPacket() noexcept {
        return Derived::isFastPacket();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval size_t device_obj<RValueVector<Derived, ScalarT>>::getSizeAtCompile() noexcept {
        return Derived::getSizeAtCompile();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval size_t device_obj<RValueVector<Derived, ScalarT>>::getSizeAtCompile(const Vector auto& hint) noexcept {
        return Derived::getSizeAtCompile(hint);
    }
    // Redeclare to expose it to the derived class
    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval void device_obj<RValueVector<Derived, ScalarT>>::static_assert_assign(const Scalar auto& source) noexcept {
        Derived::static_assert_assign(source);
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval void device_obj<RValueVector<Derived, ScalarT>>::static_assert_assign(const Vector auto& source) noexcept {
        static_assert(getSizeAtCompile() != Dynamic || DeviceObj<decltype(source)>, "[Error]: Host object cannot be assigned to dynamic device object");
        Derived::static_assert_assign(source);
    }
}
