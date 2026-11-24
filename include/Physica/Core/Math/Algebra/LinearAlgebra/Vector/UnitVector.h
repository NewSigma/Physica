/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2024 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <cstddef>
#include "Physica/Core/Scalar/Scalar.h"
#include "Vector.h"

namespace Physica {
    template<Scalar T, size_t Length = Dynamic>
    class UnitVector : public RValueVector<UnitVector<T, Length>> {
        using This = UnitVector<T, Length>;
        using Base = RValueVector<This>;

        size_t nonzero;
        size_t length;
    public:
        UnitVector(size_t nonzero_, size_t length_);
        UnitVector(const This&) = default;
        UnitVector(This&&) noexcept = default;
        ~UnitVector() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] T calc(size_t index) const noexcept;
        [[nodiscard]] constexpr size_t getNonZero() const noexcept { return nonzero; }
        [[nodiscard]] constexpr size_t getLength() const noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Length; }
    };

    template<Scalar T, size_t Length>
    UnitVector<T, Length>::UnitVector(size_t nonzero_, size_t length_) : nonzero(nonzero_), length(length_) {
        assert(nonzero < length);
    }

    template<Scalar T, size_t Length>
    void UnitVector<T, Length>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(nonzero, obj.nonzero);
        std::swap(length, obj.length);
    }

    template<Scalar T, size_t Length>
    T UnitVector<T, Length>::calc(size_t index) const noexcept {
        return T(index == nonzero ? 1 : 0);
    }

    template<Scalar T, size_t Length>
    constexpr size_t UnitVector<T, Length>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return length;
        return Length;
    }
}

namespace Physica {
    template<Scalar T, size_t Length>
    class Traits<UnitVector<T, Length>> {
    public:
        using ScalarType = T;
    };
}
