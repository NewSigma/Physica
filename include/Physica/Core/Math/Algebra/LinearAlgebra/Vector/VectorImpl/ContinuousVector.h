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

#include "Physica/Core/IO/HDF5/HDF5.h"
#include "ContinuousVectorImpl/ContinuousVectorBlock.h"

namespace Physica::Core {
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
        using typename Base::ValueType;
        using typename Base::PacketType;
        using Base::SizeAtCompile;
        using Base::isForwardDiff;
        using Base::isReverseDiff;
    protected:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    private:
        constexpr static int DiffOrder = ScalarType::Order;
        constexpr static int DataDim = 1 + (DiffOrder > 0);
        using DataSetType = H5DataSet<DataDim>;
        using DataSpaceType = H5DataSpace<DataDim>;

        using ValuesRtnTy = std::conditional<Diffable<ScalarType>, ValueVector<This>, const Derived&>::type;
    public:
        ~ContinuousVector() = default;
        /* Operators */
        using Base::operator=;
        using Base::operator+=;
        inline This& operator=(const This& v);
        inline This& operator=(This&& v);

        template<Vector V>
        inline void operator+=(const V& v);
        /* Operations */
        template<class AnyPacket> [[nodiscard]] inline AnyPacket packet(size_t index) const;
        template<class AnyPacket> [[nodiscard]] inline AnyPacket packetPartial(size_t index, size_t count) const;
        template<class AnyPacket> inline void writePacket(size_t index, const AnyPacket packet);
        template<class AnyPacket> inline void writePacketPartial(size_t index, size_t count, const AnyPacket packet);

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

        template<Vector V>
        void resize(const V& x) { resize(x.getLength()); }
        void resize(size_t length) { Base::getDerived().resize(length); }

        template<RandomGenerator R>
        inline void random_uniform();
        template<RandomGenerator R>
        inline void random_normal();
        template<class Distribution, RandomGenerator R>
        inline void random_any(Distribution& dist);

        const DataSetType read(const H5Location& loc, const char* name);
        DataSetType write(H5Location& loc, const char* name) const;

        ValuesRtnTy values() const noexcept;
        template<int GradOrder = 1>
        auto grads() const noexcept;
        /* Getters */
        [[nodiscard]] PtrTy data() { return Base::data_ptr(0); }
        [[nodiscard]] ConstPtrTy data() const { return Base::data_ptr(0); }
    protected:
        ContinuousVector() = default;
        ContinuousVector(const This&) = default;
        ContinuousVector(This&&) noexcept = default;
    };
}

#include "ContinuousVectorImpl/ContinuousVectorImpl.h"
#include "ContinuousVectorImpl/VectorConvert.h"
