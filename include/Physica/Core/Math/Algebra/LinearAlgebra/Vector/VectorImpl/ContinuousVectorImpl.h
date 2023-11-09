/*
 * Copyright 2022 WeiBo He.
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

#include "Physica/Core/Exception/NotImplementedException.h"

namespace Physica::Core {
    namespace Internal {
        template<class T1, class T2, bool enableSIMD>
        struct AddAssignImpl {
            static void run(LValueVector<T1>& v1, const RValueVector<T2>& v2) {
                using ScalarType = typename T1::ScalarType;
                for(size_t i = 0; i < v1.getLength(); ++i)
                    v1.getDerived()[i] += ScalarType(v2.calc(i));
            }
        };

        template<class T1, class T2>
        class AddAssignImpl<T1, T2, true> {
            static_assert(std::is_same<typename T1::ScalarType, typename T2::ScalarType>::value, "[Error]: SIMD on different scalars is not supported");
            static_assert(!T1::isComplex && !T2::isComplex, "[Error]: SIMD on complex scalars are not supported");
            constexpr static size_t size1 = T1::SizeAtCompile;
            constexpr static size_t size2 = T2::SizeAtCompile;
            constexpr static size_t SizeAtCompile = size1 > size2 ? size1 : size2;
            using PacketType = typename Internal::BestPacket<typename T1::ScalarType, SizeAtCompile>::Type;
        public:
            static void run(LValueVector<T1>& v1, const RValueVector<T2>& v2) {
                if constexpr (SizeAtCompile != Dynamic) {
                    constexpr size_t to = SizeAtCompile / PacketType::size() * PacketType::size();
                    for (size_t i = 0; i < to; i += PacketType::size()) {
                        const PacketType sum = v1.getDerived().template packet<PacketType>(i) + v2.getDerived().template packet<PacketType>(i);
                        v1.getDerived().writePacket(i, sum);
                    }

                    constexpr size_t i = SizeAtCompile - SizeAtCompile % PacketType::size();
                    if constexpr (i != SizeAtCompile) {
                        constexpr size_t count = SizeAtCompile - i;
                        const PacketType sum = v1.getDerived().template packetPartial<PacketType>(i, count) + v2.getDerived().template packetPartial<PacketType>(i, count);
                        v1.getDerived().writePacketPartial(i, count, sum);
                    }
                }
                else {
                    const size_t length = v1.getLength();
                    size_t i = 0;
                    const size_t to = length / PacketType::size() * PacketType::size();
                    for (; i < to; i += PacketType::size()) {
                        const PacketType sum = v1.getDerived().template packet<PacketType>(i) + v2.getDerived().template packet<PacketType>(i);
                        v1.getDerived().writePacket(i, sum);
                    }
                    if (to != length) {
                        const size_t count = length - i;
                        const PacketType sum = v1.getDerived().template packetPartial<PacketType>(i, count) + v2.getDerived().template packetPartial<PacketType>(i, count);
                        v1.getDerived().writePacketPartial(i, count, sum);
                    }
                }
            }
        };
    }

    template<class Derived>
    ContinuousVector<Derived>& ContinuousVector<Derived>::operator=(const ContinuousVector<Derived>& v) {
        Base::operator=(v);
        return *this;
    }

    template<class Derived>
    ContinuousVector<Derived>& ContinuousVector<Derived>::operator=(ContinuousVector<Derived>&& v) noexcept {
        Base::operator=(std::forward<Base>(v));
        return *this;
    }

    template<class Derived>
    template<class AnyScalar>
    inline Derived& ContinuousVector<Derived>::operator=(const ScalarBase<AnyScalar>& s) {
        static_assert(ScalarType::isComplex || !AnyScalar::isComplex, "[Error]: Assigning a complex number to real vector is not allowed");
        for (size_t i = 0; i < Base::getLength(); ++i) {
            if constexpr (AnyScalar::isReverseDiff)
                (*this)[i] = ScalarType(s.getDerived().copy());
            else
                (*this)[i] = ScalarType(s.getDerived());
        }
        return Base::getDerived();
    }

    template<class Derived>
    template<class PacketType>
    inline PacketType ContinuousVector<Derived>::packet(size_t index) const {
        assert(index + PacketType::size() <= Base::getLength());
        if constexpr (std::is_same_v<ScalarType, typename Internal::Traits<PacketType>::ScalarType>) {
            PacketType packet{};
            packet.load(Base::data_ptr(index));
            return packet;
        }
        else
            return Base::template packet<PacketType>(index);
    }

    template<class Derived>
    template<class PacketType>
    inline PacketType ContinuousVector<Derived>::packetPartial(size_t index, size_t count) const {
        if constexpr (std::is_same_v<ScalarType, typename Internal::Traits<PacketType>::ScalarType>) {
            PacketType packet{};
            packet.load_partial(count, Base::data_ptr(index));
            return packet;
        }
        else
            return Base::template packetPartial<PacketType>(index, count);
    }

    template<class Derived>
    template<class PacketType>
    inline void ContinuousVector<Derived>::writePacket(size_t index, const PacketType packet) {
        constexpr bool isSameScalar = std::is_same_v<ScalarType, typename Internal::Traits<PacketType>::ScalarType>;
        if constexpr (isSameScalar)
            packet.store(Base::data_ptr(index));
        else
            Base::template writePacket(index, packet);
    }

    template<class Derived>
    template<class PacketType>
    inline void ContinuousVector<Derived>::writePacketPartial(size_t index, size_t count, const PacketType packet) {
        constexpr bool isSameScalar = std::is_same_v<ScalarType, typename Internal::Traits<PacketType>::ScalarType>;
        if constexpr (isSameScalar)
            packet.store_partial(count, Base::data_ptr(index));
        else
            Base::template writePacketPartial<PacketType>(index, count, packet);
    }

    template<class Derived>
    template<size_t Length>
    inline ContinuousVectorBlock<Derived, Length> ContinuousVector<Derived>::head(size_t to) {
        return ContinuousVectorBlock<Derived, Length>(Base::getDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    inline const ContinuousVectorBlock<Derived, Length> ContinuousVector<Derived>::head(size_t to) const {
        return ContinuousVectorBlock<Derived, Length>(Base::getConstCastDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    inline ContinuousVectorBlock<Derived, Length> ContinuousVector<Derived>::tail(size_t from) {
        return ContinuousVectorBlock<Derived, Length>(Base::getDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    inline const ContinuousVectorBlock<Derived, Length> ContinuousVector<Derived>::tail(size_t from) const {
        return ContinuousVectorBlock<Derived, Length>(Base::getConstCastDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    inline ContinuousVectorBlock<Derived, Length> ContinuousVector<Derived>::segment(size_t from, size_t to) {
        return ContinuousVectorBlock<Derived, Length>(Base::getDerived(), from, to);
    }

    template<class Derived>
    template<size_t Length>
    inline const ContinuousVectorBlock<Derived, Length> ContinuousVector<Derived>::segment(size_t from, size_t to) const {
        return ContinuousVectorBlock<Derived, Length>(Base::getConstCastDerived(), from, to);
    }

    template<class Derived>
    bool ContinuousVector<Derived>::checkContinuous() const {
        if constexpr (isReverseDiff) {
            const size_t traceIndex0 = (*this)[0].getTraceIndex();
            for (size_t i = 1; i < Base::getLength(); ++i)
                if ((*this)[i].getTraceIndex() != traceIndex0 + i)
                    return false;
        }
        return true;
    }

    template<class Derived>
    void ContinuousVector<Derived>::makeContinuous() {
        if constexpr (isReverseDiff) {
            for (size_t i = 0; i < Base::getLength(); ++i) {
                auto& elem = (*this)[i];
                elem = elem.copy();
            }
        }
    }

    template<class Derived>
    template<class SpaceType>
    void ContinuousVector<Derived>::read(const H5::DataSet& dataset,
                                         const DataSpaceBase<SpaceType>& file_space,
                                         const H5::DSetMemXferPropList& xfer_plist) {
        const auto mem_space = H5DataSpace<1>({Base::getLength()});
        dataset.read(data(), ScalarType::getH5DataType(), mem_space, file_space.asH5Type(), xfer_plist);
    }

    template<class Derived>
    template<class SpaceType>
    void ContinuousVector<Derived>::write(H5::DataSet& dataset,
                                          const DataSpaceBase<SpaceType>& file_space,
                                          const H5::DSetMemXferPropList& xfer_plist) const {
        const auto mem_space = H5DataSpace<1>({Base::getLength()});
        dataset.write(data(), ScalarType::getH5DataType(), mem_space, file_space.asH5Type(), xfer_plist);
    }

    template<class Derived>
    std::ostream& operator<<(std::ostream& os, const ContinuousVector<Derived>& v) {
        using ScalarType = typename Derived::ScalarType;
        os.write(reinterpret_cast<const char*>(v.data_ptr(0)), v.getLength() * sizeof(ScalarType));
        return os;
    }

    template<class Derived>
    std::istream& operator>>(std::istream& is, ContinuousVector<Derived>& v) {
        using ScalarType = typename Derived::ScalarType;
        is.read(reinterpret_cast<char*>(v.data_ptr(0)), v.getLength() * sizeof(ScalarType));
        return is;
    }

    template<class Derived, class OtherDerived>
    void operator+=(ContinuousVector<Derived>& v1, const RValueVector<OtherDerived>& v2) {
        assert(v1.getLength() == v2.getLength());
        Internal::AddAssignImpl<Derived, OtherDerived, Internal::EnableSIMD<Derived, OtherDerived>::value>::run(v1, v2);
    }
}
