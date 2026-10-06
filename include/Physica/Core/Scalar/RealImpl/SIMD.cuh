/*
 * Copyright 2026 Weibo He.
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

#include "Physica/Core/Scalar/Real.h"
#include "Physica/Core/Scalar/ScalarImpl/SIMDMixin.h"
#include "Physica/Core/Utils/CUDA/device_obj.h"
#include "Physica/Core/Utils/Empty.h"

namespace Physica {
    template<Scalar T, int Size>
    class device_obj<SIMD<T, Size>> : public SIMDMixin<device_obj<SIMD<T, Size>>> {
        using This = device_obj<SIMD<T, Size>>;
        using Base = SIMDMixin<This>;
        using Pack = Traits<This>::MachineType;
        constexpr static bool isHalf2 = std::same_as<Pack, __half2>;
        constexpr static bool isFloat4 = std::same_as<Pack, float4>;

        static_assert(Either<Pack, __half2, float4, double2>, "[Error]: Unexpected pack type");
    private:
        Pack pack;
    public:
        constexpr device_obj() noexcept = default;
        __device__ explicit device_obj(T x) noexcept;
        __device__ device_obj(T x, int count) noexcept;
        __device__ explicit device_obj(Pack x) noexcept;
        constexpr device_obj(const This&) = default;
        constexpr device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        [[nodiscard]] __device__ T operator[](int index) const noexcept;
        [[nodiscard]] __device__ This operator+(const This& other) const noexcept;
        [[nodiscard]] __device__ This operator-(const This& other) const noexcept;
        [[nodiscard]] __device__ This operator*(const This& other) const noexcept;
        [[nodiscard]] __device__ This operator*(const T& x) const noexcept;
        [[nodiscard]] __device__ This operator/(const This& other) const noexcept;
        [[nodiscard]] __device__ This operator-() const noexcept;
        __device__ void operator+=(const This& other) noexcept;
        __device__ void operator-=(const This& other) noexcept;
        __device__ void operator*=(const This& other) noexcept;
        __device__ void operator*=(const T& x) noexcept;
        __device__ void operator/=(const This& other) noexcept;
        /* Operations */
        __device__ void load(const T* p) noexcept;
        __device__ void load(const T* p, int n) noexcept;
        __device__ void store(T* p) const noexcept;
        __device__ void store(T* p, int n) const noexcept;
        __device__ void insert(int index, const T& value) noexcept;

        [[nodiscard]] __device__ T sum() const noexcept;
        [[nodiscard]] __device__ T max() const noexcept;
        [[nodiscard]] __device__ T min() const noexcept;
        /* Getters */
        [[nodiscard]] __device__ Pack& toMachine() noexcept;
        [[nodiscard]] __device__ const Pack& toMachine() const noexcept;
    };

    template<Scalar T, int Size>
    __device__ device_obj<SIMD<T, Size>>::device_obj(const T x) noexcept {
        const auto lane = x.toMachine();
        if constexpr (isHalf2)
            pack = __halves2half2(lane, lane);
        else if constexpr (isFloat4)
            pack = make_float4(lane, lane, lane, lane);
        else
            pack = make_double2(lane, lane);
    }

    template<Scalar T, int Size>
    __device__ device_obj<SIMD<T, Size>>::device_obj(const T x, const int count) noexcept {
        assert(0 < count && count <= int(Size) && "[Error]: Invalid count");
        const auto lane = x.toMachine();
        const auto zero = T(0).toMachine();
        if constexpr (isHalf2)
            pack = __halves2half2(lane, count > 1 ? lane : zero);
        else if constexpr (isFloat4) {
            switch (count) {
            case 1:
                pack = make_float4(lane, zero, zero, zero);
                break;
            case 2:
                pack = make_float4(lane, lane, zero, zero);
                break;
            case 3:
                pack = make_float4(lane, lane, lane, zero);
                break;
            case 4:
                pack = make_float4(lane, lane, lane, lane);
                break;
            default:
                unreachable();
            }
        }
        else
            pack = make_double2(lane, count > 1 ? lane : zero);
    }

    template<Scalar T, int Size>
    __device__ device_obj<SIMD<T, Size>>::device_obj(const Pack x) noexcept : pack(x) {}

    template<Scalar T, int Size>
    __device__ auto device_obj<SIMD<T, Size>>::operator[](const int index) const noexcept -> T {
        assert(0 <= index && index < Size && "[Error]: Index out of range");
        if constexpr (isHalf2)
            return T(index == 0 ? __low2half(pack) : __high2half(pack));
        else if constexpr (isFloat4) {
            switch (index) {
            case 0:
                return pack.x;
            case 1:
                return pack.y;
            case 2:
                return pack.z;
            case 3:
                return pack.w;
            default:
                unreachable();
            }
        }
        else
            return T(index == 0 ? pack.x : pack.y);
    }

    template<Scalar T, int Size>
    __device__ auto device_obj<SIMD<T, Size>>::operator+(const This& other) const noexcept -> This {
        if constexpr (isHalf2)
            return This(__hadd2(pack, other.pack));
        else if constexpr (isFloat4)
            return This(make_float4(pack.x + other.pack.x, pack.y + other.pack.y, pack.z + other.pack.z, pack.w + other.pack.w));
        else
            return This(make_double2(pack.x + other.pack.x, pack.y + other.pack.y));
    }

    template<Scalar T, int Size>
    __device__ auto device_obj<SIMD<T, Size>>::operator-(const This& other) const noexcept -> This {
        if constexpr (isHalf2)
            return This(__hsub2(pack, other.pack));
        else if constexpr (isFloat4)
            return This(make_float4(pack.x - other.pack.x, pack.y - other.pack.y, pack.z - other.pack.z, pack.w - other.pack.w));
        else
            return This(make_double2(pack.x - other.pack.x, pack.y - other.pack.y));
    }

    template<Scalar T, int Size>
    __device__ auto device_obj<SIMD<T, Size>>::operator*(const This& other) const noexcept -> This {
        if constexpr (isHalf2)
            return This(__hmul2(pack, other.pack));
        else if constexpr (isFloat4)
            return This(make_float4(pack.x * other.pack.x, pack.y * other.pack.y, pack.z * other.pack.z, pack.w * other.pack.w));
        else
            return This(make_double2(pack.x * other.pack.x, pack.y * other.pack.y));
    }

    template<Scalar T, int Size>
    __device__ auto device_obj<SIMD<T, Size>>::operator*(const T& x) const noexcept -> This {
        return operator*(This(x));
    }

    template<Scalar T, int Size>
    __device__ auto device_obj<SIMD<T, Size>>::operator/(const This& other) const noexcept -> This {
        if constexpr (isHalf2)
            return This(__h2div(pack, other.pack));
        else if constexpr (isFloat4)
            return This(make_float4(pack.x / other.pack.x, pack.y / other.pack.y, pack.z / other.pack.z, pack.w / other.pack.w));
        else
            return This(make_double2(pack.x / other.pack.x, pack.y / other.pack.y));
    }

    template<Scalar T, int Size>
    __device__ auto device_obj<SIMD<T, Size>>::operator-() const noexcept -> This {
        if constexpr (isHalf2)
            return This(__hneg2(pack));
        else if constexpr (isFloat4)
            return This(make_float4(-pack.x, -pack.y, -pack.z, -pack.w));
        else
            return This(make_double2(-pack.x, -pack.y));
    }

    template<Scalar T, int Size>
    __device__ void device_obj<SIMD<T, Size>>::operator+=(const This& other) noexcept {
        *this = *this + other;
    }

    template<Scalar T, int Size>
    __device__ void device_obj<SIMD<T, Size>>::operator-=(const This& other) noexcept {
        *this = *this - other;
    }

    template<Scalar T, int Size>
    __device__ void device_obj<SIMD<T, Size>>::operator*=(const This& other) noexcept {
        *this = *this * other;
    }

    template<Scalar T, int Size>
    __device__ void device_obj<SIMD<T, Size>>::operator*=(const T& x) noexcept {
        *this = *this * x;
    }

    template<Scalar T, int Size>
    __device__ void device_obj<SIMD<T, Size>>::operator/=(const This& other) noexcept {
        *this = *this / other;
    }

    template<Scalar T, int Size>
    __device__ void device_obj<SIMD<T, Size>>::load(const T* p) noexcept {
        assert(p % alignof(Pack) == 0 && "[Error]: Unaligned load");
        pack = *reinterpret_cast<const Pack*>(p);
    }

    template<Scalar T, int Size>
    __device__ void device_obj<SIMD<T, Size>>::load(const T* p, const int n) noexcept {
        assert(0 < n && n <= int(Size) && "[Error]: Invalid count for partial load");
        *this = This(T(0));
        for (int i = 0; i < n; ++i)
            insert(i, p[i]);
    }

    template<Scalar T, int Size>
    __device__ void device_obj<SIMD<T, Size>>::store(T* p) const noexcept {
        assert(p % alignof(Pack) == 0 && "[Error]: Unaligned store");
        *reinterpret_cast<Pack*>(p) = pack;
    }

    template<Scalar T, int Size>
    __device__ void device_obj<SIMD<T, Size>>::store(T* p, const int n) const noexcept {
        assert(0 < n && n <= int(Size) && "[Error]: Invalid count for partial store");
        for (int i = 0; i < n; ++i)
            p[i] = (*this)[i];
    }

    template<Scalar T, int Size>
    __device__ void device_obj<SIMD<T, Size>>::insert(const int index, const T& value) noexcept {
        assert(0 <= index && index < Size && "[Error]: Index out of range");
        const auto lane = value.toMachine();
        if constexpr (isHalf2) {
            switch (index) {
            case 0:
                pack = __halves2half2(lane, __high2half(pack));
                break;
            case 1:
                pack = __halves2half2(__low2half(pack), lane);
                break;
            default:
                unreachable();
            }
        }
        else if constexpr (isFloat4) {
            switch (index) {
            case 0:
                pack = make_float4(lane, pack.y, pack.z, pack.w);
                break;
            case 1:
                pack = make_float4(pack.x, lane, pack.z, pack.w);
                break;
            case 2:
                pack = make_float4(pack.x, pack.y, lane, pack.w);
                break;
            case 3:
                pack = make_float4(pack.x, pack.y, pack.z, lane);
                break;
            default:
                unreachable();
            }
        }
        else {
            switch (index) {
            case 0:
                pack = make_double2(lane, pack.y);
                break;
            case 1:
                pack = make_double2(pack.x, lane);
                break;
            default:
                unreachable();
            }
        }
    }

    template<Scalar T, int Size>
    __device__ auto device_obj<SIMD<T, Size>>::sum() const noexcept -> T {
        T result{};
        for (int i = 0; i < Size; ++i)
            result = result + (*this)[i];
        return result;
    }

    template<Scalar T, int Size>
    __device__ auto device_obj<SIMD<T, Size>>::max() const noexcept -> T {
        T result = (*this)[0];
        for (int i = 1; i < Size; ++i)
            result = std::max(result, (*this)[i]);
        return result;
    }

    template<Scalar T, int Size>
    __device__ auto device_obj<SIMD<T, Size>>::min() const noexcept -> T {
        T result = (*this)[0];
        for (int i = 1; i < Size; ++i)
            result = std::min(result, (*this)[i]);
        return result;
    }

    template<Scalar T, int Size>
    __device__ auto device_obj<SIMD<T, Size>>::toMachine() noexcept -> Pack& {
        return pack;
    }

    template<Scalar T, int Size>
    __device__ auto device_obj<SIMD<T, Size>>::toMachine() const noexcept -> const Pack& {
        return pack;
    }
}

namespace Physica {
    template<Scalar T, int N>
    class Traits<device_obj<SIMD<T, N>>> {
        template<class> struct PacketType;

        template<>
        struct PacketType<Real<Float16>> {
            constexpr static int Size = 2;
            using Type = __half2;
        };

        template<>
        struct PacketType<Real<Float32>> {
            constexpr static int Size = 4;
            using Type = float4;
        };

        template<>
        struct PacketType<Real<Float64>> {
            constexpr static int Size = 2;
            using Type = double2;
        };
    public:
        constexpr static int Size = PacketType<T>::Size;
        static_assert(N == Size, "[Error]: Device packet must be the natural width of its scalar");
        using ScalarType = T;
        using ValueType = device_obj<SIMD<T, N>>;
        using GradType = void;
        using RealType = ValueType;
        using FullRealType = ValueType;
        using MachineType = PacketType<T>::Type;
        using BoolSIMDType = Empty;
        constexpr static bool isSeparatable = false;
    };
}

#include "SIMDImpl/BestPacket.cuh"
#include "SIMDImpl/Math.cuh"
