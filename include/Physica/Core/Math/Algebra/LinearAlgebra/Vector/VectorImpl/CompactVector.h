/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "CompactVectorImpl/CompactVectorBlock.h"
#include "StridedVector.h"

namespace Physica {
    template<class Derived>
    class CompactVector : public StridedVector<Derived> {
        using Base = StridedVector<Derived>;
        using This = CompactVector<Derived>;
    public:
        using Base::isDiffable;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        constexpr static int DiffOrder = T::Order;
        constexpr static int DataDim = 1 + (DiffOrder > 0);
        using DataSetType = H5Dataset<DataDim>;
        using DataSpaceType = H5DataSpace<DataDim>;

        template<Vector> class View;
    public:
        ~CompactVector() = default;
        /* Operators */
        This& operator=(const This& v) = delete;
        This& operator=(This&& v) noexcept = delete;
        [[nodiscard]] bool operator==(this const auto& self, const Vector auto& other) noexcept;
        using Base::operator=;
        /* Operations */
        template<int Size> [[nodiscard]] SIMD<T, Size> packet(size_t index) const noexcept;
        template<int Size> [[nodiscard]] SIMD<T, Size> packet(size_t index, size_t count) const noexcept;
        void writePacket(Packet auto packet, size_t index) noexcept;
        void writePacket(Packet auto packet, size_t index, size_t count) noexcept;
        [[nodiscard]] constexpr auto view(this auto&&) noexcept;

        template<Vector V> void toDevice(device_obj<CompactVector<V>>& obj) const;
        template<Vector V> void toDeviceAsync(device_obj<CompactVector<V>>& obj) const;

        template<size_t Length = Dynamic>
        [[nodiscard]] auto head(this auto&&, size_t to = Length) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] auto tail(this auto&&, size_t from) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] auto segment(this auto&&, size_t from, size_t to) noexcept;

        void read(const auto& obj);
        const DataSetType read(const H5Loc& loc, const char* name);
        DataSetType write(H5Loc& loc, const char* name) const;
        /* Getters */
        [[nodiscard]] constexpr auto data(this auto&&) noexcept;
        [[nodiscard]] auto data_handle() noexcept;
        [[nodiscard]] auto data_handle() const noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getStrideAtCompile() noexcept { return 1; }
    protected:
        CompactVector() = default;
        CompactVector(const This&) = default;
        CompactVector(This&&) noexcept = default;
        /* Friends */
        friend class device_obj<CompactVector<Derived>>;
    };
}

#include "CompactVectorImpl/CompactVectorImpl.h"
#include "CompactVectorImpl/View.h"
