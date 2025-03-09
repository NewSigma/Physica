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

#include "ContinuousVectorImpl/ContinuousVectorBlock.h"

namespace Physica {
    /**
     * \class ContinuousVector has its element on memory continuously
     */
    template<class Derived>
    class ContinuousVector : public LValueVector<Derived> {
        using Base = LValueVector<Derived>;
        using This = ContinuousVector<Derived>;
        template<size_t Length>
        using BlockType = ContinuousVectorBlock<Derived, Length>;
    public:
        using typename Base::ScalarType;
        using typename Base::PacketType;
        using Base::SizeAtCompile;
        using Base::isForwardDiff;
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    private:
        constexpr static int DiffOrder = ScalarType::Order;
        constexpr static int DataDim = 1 + (DiffOrder > 0);
        using DataSetType = H5DataSet<DataDim>;
        using DataSpaceType = H5DataSpace<DataDim>;

        using ValuesRtnTy = std::conditional<Diffable<ScalarType>, ValueVector<This>, Derived&>::type;
    public:
        ~ContinuousVector() = default;
        /* Operators */
        using Base::operator=;
        using Base::operator+=;
        inline This& operator=(const This& v);
        inline This& operator=(This&& v);
        /* Operations */
        template<Packet Pack> [[nodiscard]] inline Pack packet(size_t index) const;
        template<Packet Pack> [[nodiscard]] inline Pack packetPartial(size_t index, size_t count) const;
        template<Packet Pack> inline void writePacket(size_t index, const Pack packet);
        template<Packet Pack> inline void writePacketPartial(size_t index, size_t count, const Pack packet);

        template<class U>
        void reverse(const U& grad) const noexcept requires(isReverseDiff);

        template<Vector V>
        void resize(const V& x) { resize(x.getLength()); }
        auto resize(size_t length) { return Base::getDerived().resize(length); }
        template<Vector T> void toDevice(device_obj<ContinuousVector<T>>& obj) const;
        template<Vector T> void toDeviceAsync(device_obj<ContinuousVector<T>>& obj) const;
        [[nodiscard]] auto toNumpy() const;

        template<size_t Length = Dynamic>
        [[nodiscard]] inline auto head(size_t to) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline const auto head(size_t to) const noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline auto tail(size_t from) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline const auto tail(size_t from) const noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline auto segment(size_t from, size_t to) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline const auto segment(size_t from, size_t to) const noexcept;

        void zeros();
        template<RNG R>
        inline void random_uniform();
        template<RNG R>
        inline void random_normal();
        template<RNG R, class Distribution>
        inline void random_any(Distribution& dist);

        const DataSetType read(const H5Loc& loc, const char* name);
        DataSetType write(H5Loc& loc, const char* name) const;

        [[nodiscard]] ValuesRtnTy values() noexcept;
        [[nodiscard]] const ValuesRtnTy values() const noexcept;
        template<int GradOrder = 1>
        auto grads() const noexcept;
        /* Getters */
        [[nodiscard]] PtrTy data() { return Base::data_ptr(0); }
        [[nodiscard]] ConstPtrTy data() const { return Base::data_ptr(0); }
    protected:
        ContinuousVector() = default;
        ContinuousVector(const This&) = default;
        ContinuousVector(This&&) noexcept = default;
        /* Friends */
        friend class device_obj<ContinuousVector<Derived>>;
    };
}

#include "ContinuousVectorImpl/ContinuousVectorImpl.h"
#include "ContinuousVectorImpl/VectorConvert.h"
