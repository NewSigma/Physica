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

#include "CompactVectorImpl/CompactVectorBlock.h"
#include "StridedVector.h"

namespace Physica {
    template<Vector, int MatrixMajor, size_t Row, size_t Col>
    class CompactReshapedVector;

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
        using DataSetType = H5DataSet<DataDim>;
        using DataSpaceType = H5DataSpace<DataDim>;

        template<Vector> class View;
    public:
        ~CompactVector() = default;
        /* Operators */
        This& operator=(const This& v) = delete;
        This& operator=(This&& v) noexcept = delete;
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
        template<int Major, size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto reshape(this auto&& self, size_t row, size_t col) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto reshape_row(this auto&& self, size_t row, size_t col) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto reshape_col(this auto&& self, size_t row, size_t col) noexcept;
        [[nodiscard]] auto reshape_like(this auto&& self, const Matrix auto& mat) noexcept;

        void read(const auto& obj) noexcept;
        const DataSetType read(const H5Loc& loc, const char* name);
        DataSetType write(H5Loc& loc, const char* name) const;
        /* Getters */
        [[nodiscard]] auto data() noexcept;
        [[nodiscard]] auto data() const noexcept;
        [[nodiscard]] auto data_handle(this auto&&) noexcept;
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
#include "CompactVectorImpl/VectorConvert/RealVector.h"
#include "CompactVectorImpl/VectorConvert/ImagVector.h"
