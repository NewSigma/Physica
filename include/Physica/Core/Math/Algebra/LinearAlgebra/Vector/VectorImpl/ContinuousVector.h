/*
 * Copyright 2022-2023 WeiBo He.
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

#include "LValueVector.h"
#include "ContinuousVectorBlock.h"
#include "Physica/Core/IO/HDF5/HDF5.h"

namespace Physica::Core {
    /**
     * \class ContinuousVector has its element on memory continuously
     */
    template<class Derived>
    class ContinuousVector : public LValueVector<Derived> {
    public:
        using Base = LValueVector<Derived>;
        using typename Base::ScalarType;
        using typename Base::PlainScalar;
        using Base::isReverseDiff;
    public:
        ~ContinuousVector() = default;
        /* Operators */
        using Base::operator=;
        inline ContinuousVector& operator=(const ContinuousVector& v);
        inline ContinuousVector& operator=(ContinuousVector&& v) noexcept;
        inline Derived& operator=(const ScalarType& s);
        /* Operations */
        template<class PacketType> [[nodiscard]] inline PacketType packet(size_t index) const;
        template<class PacketType> [[nodiscard]] inline PacketType packetPartial(size_t index, size_t count) const;
        template<class PacketType> inline void writePacket(size_t index, const PacketType packet);
        template<class PacketType> inline void writePacketPartial(size_t index, size_t count, const PacketType packet);
        void resize(size_t length) { Base::getDerived().resize(length); }

        template<class OtherDerived> void toDevice(device_obj<ContinuousVector<OtherDerived>>& obj) const;
        template<class OtherDerived> void toDeviceAsync(device_obj<ContinuousVector<OtherDerived>>& obj) const;

        template<size_t Length = Dynamic> inline ContinuousVectorBlock<Derived, Length> head(size_t to);
        template<size_t Length = Dynamic> inline const ContinuousVectorBlock<Derived, Length> head(size_t to) const;
        template<size_t Length = Dynamic> inline ContinuousVectorBlock<Derived, Length> tail(size_t from);
        template<size_t Length = Dynamic> inline const ContinuousVectorBlock<Derived, Length> tail(size_t from) const;
        template<size_t Length = Dynamic> inline ContinuousVectorBlock<Derived, Length> segment(size_t from, size_t to);
        template<size_t Length = Dynamic> inline const ContinuousVectorBlock<Derived, Length> segment(size_t from, size_t to) const;

        [[nodiscard]] bool checkContinuous() const;
        inline void makeContinuous();
        template<class RandomGenerator>
        inline void random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        inline void random_normal(RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        inline void random_any(Distribution& dist, RandomGenerator& gen);

        void read(const H5File& h5f, const char* name, const H5::DSetMemXferPropList& xfer_plist = H5::DSetMemXferPropList::DEFAULT);
        template<class SpaceType>
        void read(const H5::DataSet& dataset,
                  const DataSpaceBase<SpaceType>& file_space,
                  const H5::DSetMemXferPropList& xfer_plist = H5::DSetMemXferPropList::DEFAULT);
        
        void write(H5File& h5f, const char* name, const H5::DSetMemXferPropList& xfer_plist = H5::DSetMemXferPropList::DEFAULT) const;
        template<class SpaceType>
        void write(H5::DataSet& dataset, 
                   const DataSpaceBase<SpaceType>& file_space,
                   const H5::DSetMemXferPropList& xfer_plist = H5::DSetMemXferPropList::DEFAULT) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ ScalarType* data() { return Base::data_ptr(0); }
        [[nodiscard]] __host__ __device__ const ScalarType* data() const { return Base::data_ptr(0); }
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
    void operator+=(ContinuousVector<Derived>& v1, const RValueVector<OtherDerived>& v2);
}

#include "ContinuousVectorImpl.h"
