/*
 * Copyright 2022-2024 Weibo He.
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

#include <Physica/Core/IO/HDF5/HDF5.h>
#include "LValueVector.h"
#include "ContinuousVectorImpl/ContinuousVectorBlock.h"

namespace Physica::Core {
    /**
     * \class ContinuousVector has its element on memory continuously
     */
    template<class Derived>
    class ContinuousVector : public LValueVector<Derived> {
        using Base = LValueVector<Derived>;
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
        constexpr static int DataDim = 1 + (DiffOrder > 1);
        using DataSetType = H5DataSet<DataDim>;
        using DataSpaceType = H5DataSpace<DataDim>;
    public:
        ~ContinuousVector() = default;
        /* Operators */
        using Base::operator=;
        inline ContinuousVector& operator=(const ContinuousVector& v);
        inline ContinuousVector& operator=(ContinuousVector&& v) noexcept;
        template<Scalar T>
        inline Derived& operator=(const T& x);
        /* Operations */
        template<class AnyPacket> [[nodiscard]] inline AnyPacket packet(size_t index) const;
        template<class AnyPacket> [[nodiscard]] inline AnyPacket packetPartial(size_t index, size_t count) const;
        template<class AnyPacket> inline void writePacket(size_t index, const AnyPacket packet);
        template<class AnyPacket> inline void writePacketPartial(size_t index, size_t count, const AnyPacket packet);

        template<class OtherDerived> void toDevice(device_obj<ContinuousVector<OtherDerived>>& obj) const;
        template<class OtherDerived> void toDeviceAsync(device_obj<ContinuousVector<OtherDerived>>& obj) const;

        template<size_t Length = Dynamic>
        [[nodiscard]] inline BlockType<Length> head(size_t to);
        template<size_t Length = Dynamic>
        [[nodiscard]] inline const BlockType<Length> head(size_t to) const;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline BlockType<Length> tail(size_t from);
        template<size_t Length = Dynamic>
        [[nodiscard]] inline const BlockType<Length> tail(size_t from) const;
        template<size_t Length = Dynamic>
        [[nodiscard]] inline BlockType<Length> segment(size_t from, size_t to);
        template<size_t Length = Dynamic>
        [[nodiscard]] inline const BlockType<Length> segment(size_t from, size_t to) const;

        void resize(size_t length) { Base::getDerived().resize(length); }
        inline void makeContinuous();
        [[nodiscard]] bool checkContinuous() const;
        template<class RandomGenerator>
        inline void random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        inline void random_normal(RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        inline void random_any(Distribution& dist, RandomGenerator& gen);

        const DataSetType read(const H5Location& loc, const char* name);
        DataSetType write(H5Location& loc, const char* name) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ PtrTy data() { return Base::data_ptr(0); }
        [[nodiscard]] __host__ __device__ ConstPtrTy data() const { return Base::data_ptr(0); }
    protected:
        ContinuousVector() = default;
        ContinuousVector(const ContinuousVector&) = default;
        ContinuousVector(ContinuousVector&&) noexcept = default;
    private:
        /* Friends */
        template<class T>
        friend std::ostream& operator<<(std::ostream& os, const ContinuousVector<T>& v);
        template<class T>
        friend std::istream& operator>>(std::istream& is, ContinuousVector<T>& v);
    };

    template<class Derived>
    std::ostream& operator<<(std::ostream& os, const ContinuousVector<Derived>& v);

    template<class Derived>
    std::istream& operator>>(std::istream& is, ContinuousVector<Derived>& v);

    template<class Derived, class OtherDerived>
    inline void operator+=(ContinuousVector<Derived>& v1, const RValueVector<OtherDerived>& v2);
}

#include "ContinuousVectorImpl/ContinuousVectorImpl.h"
