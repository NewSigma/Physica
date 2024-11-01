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

#include <Physica/Core/Exception/NoImplException.h>
#include <Physica/Core/Exception/MKL/VSL.h>

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
            constexpr static size_t size1 = T1::SizeAtCompile;
            constexpr static size_t size2 = T2::SizeAtCompile;
            constexpr static size_t SizeAtCompile = size1 > size2 ? size1 : size2;
            using PacketType = typename BestPacket<typename T1::ScalarType, SizeAtCompile>::Type;
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
    inline ContinuousVector<Derived>& ContinuousVector<Derived>::operator=(const ContinuousVector<Derived>& v) {
        Base::operator=(v);
        return *this;
    }

    template<class Derived>
    inline ContinuousVector<Derived>& ContinuousVector<Derived>::operator=(ContinuousVector<Derived>&& v) noexcept {
        Base::operator=(std::forward<Base>(v));
        return *this;
    }

    template<class Derived>
    inline Derived& ContinuousVector<Derived>::operator=(const ScalarType& s) {
        const size_t length = Base::getLength();
        if constexpr (isReverseDiff) {
            using TracerType = typename ScalarType::TracerType;
            TracerType::getInstance().reserve(length);
            for (size_t i = 0; i < length; ++i)
                (*this)[i] = s.copy();
        }
        else {
            for (size_t i = 0; i < length; ++i)
                (*this)[i] = s;
        }
        return Base::getDerived();
    }

    template<class Derived>
    template<class AnyPacket>
    inline AnyPacket ContinuousVector<Derived>::packet(size_t index) const {
        assert(index + AnyPacket::size() <= Base::getLength());
        if constexpr (std::is_same_v<ScalarType, typename Traits<AnyPacket>::ScalarType>) {
            AnyPacket packet{};
            packet.load(Base::data_ptr(index));
            return packet;
        }
        else
            return Base::template packet<AnyPacket>(index);
    }

    template<class Derived>
    template<class AnyPacket>
    inline AnyPacket ContinuousVector<Derived>::packetPartial(size_t index, size_t count) const {
        assert(index + count <= Base::getLength());
        if constexpr (std::is_same_v<ScalarType, typename Traits<AnyPacket>::ScalarType>) {
            AnyPacket packet{};
            packet.load_partial(count, Base::data_ptr(index));
            return packet;
        }
        else
            return Base::template packetPartial<AnyPacket>(index, count);
    }

    template<class Derived>
    template<class AnyPacket>
    inline void ContinuousVector<Derived>::writePacket(size_t index, const AnyPacket packet) {
        constexpr bool isSameScalar = std::is_same_v<ScalarType, typename Traits<AnyPacket>::ScalarType>;
        if constexpr (isSameScalar)
            packet.store(Base::data_ptr(index));
        else
            Base::template writePacket(index, packet);
    }

    template<class Derived>
    template<class AnyPacket>
    inline void ContinuousVector<Derived>::writePacketPartial(size_t index, size_t count, const AnyPacket packet) {
        constexpr bool isSameScalar = std::is_same_v<ScalarType, typename Traits<AnyPacket>::ScalarType>;
        if constexpr (isSameScalar)
            packet.store_partial(count, Base::data_ptr(index));
        else
            Base::template writePacketPartial<AnyPacket>(index, count, packet);
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
    inline void ContinuousVector<Derived>::makeContinuous() {
        if constexpr (isReverseDiff)
            Base::getDerived() = Base::getDerived().copy();
    }

    template<class Derived>
    bool ContinuousVector<Derived>::checkContinuous() const {
        if constexpr (isReverseDiff) {
            const ScalarType i0 = (*this)[0];
            const auto* pValue = i0.value_ptr();
            const auto* pGrad = i0.grad_ptr();
            for (size_t i = 1; i < Base::getLength(); ++i) {
                const ScalarType& elem = (*this)[i];
                if (elem.value_ptr() != pValue + i)
                    return false;
                if (elem.grad_ptr() != pGrad + i)
                    return false;
            }
        }
        return true;
    }

    template<class Derived>
    template<class RandomGenerator>
    inline void ContinuousVector<Derived>::random_uniform(RandomGenerator& gen) {
        if constexpr (isReverseDiff) {
            using TracerType = typename ScalarType::TracerType;
            const size_t length = Base::getLength();
            TracerType::getInstance().reserve(length);
            for (size_t i = 0; i < length; ++i)
                this->operator[](i) = ScalarType::random_uniform(gen);
        }
        else if constexpr (HasMKL()) {
            const size_t length = Base::getLength() * (Base::isComplex ? 2 : 1) * (Base::isForwardDiff ? 2 : 1);
            if constexpr (ScalarType::Option == Float32)
                vslCheck(vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, gen, length, (float*)data(), 0, 1));
            else if constexpr (ScalarType::Option == Float64)
                vslCheck(vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, gen, length, (double*)data(), 0, 1));
            else
                Base::random_uniform(gen);
        }
        else
            Base::random_uniform(gen);
    }

    template<class Derived>
    template<class RandomGenerator>
    inline void ContinuousVector<Derived>::random_normal(RandomGenerator& gen) {
        if constexpr (isReverseDiff) {
            using TracerType = typename ScalarType::TracerType;
            const size_t length = Base::getLength();
            TracerType::getInstance().reserve(length);
            for (size_t i = 0; i < length; ++i)
                this->operator[](i) = ScalarType::random_normal(gen);
        }
        else if constexpr (HasMKL() && !isForwardDiff) {
            const size_t length = Base::getLength() * (Base::isComplex ? 2 : 1);
            if constexpr (ScalarType::Option == Float32)
                vslCheck(vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, gen, length, (float*)data(), 0, 1));
            else if constexpr (ScalarType::Option == Float64)
                vslCheck(vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, gen, length, (double*)data(), 0, 1));
            else
                Base::random_normal(gen);
        }
        else
            Base::random_normal(gen);
    }

    template<class Derived>
    template<class Distribution, class RandomGenerator>
    inline void ContinuousVector<Derived>::random_any(Distribution& dist, RandomGenerator& gen) {
        if constexpr (isReverseDiff) {
            using TracerType = typename ScalarType::TracerType;
            const size_t length = Base::getLength();
            TracerType::getInstance().reserve(length);
            for (size_t i = 0; i < length; ++i)
                this->operator[](i) = ScalarType::random_any(dist, gen);
        }
        else
            Base::random_any(dist, gen);
    }
#ifdef PHYSICA_HDF5
    template<class Derived>
    const H5DataSet<1> ContinuousVector<Derived>::read(const H5Location& loc, const char* name) {
        const auto dataset = loc.openDataSet<1>(name);
        const size_t length = dataset.getSize(0);
        resize(length);
        const auto space = H5DataSpace<1>({length});
        dataset.read(data(), ScalarType::getH5DataType(), space, space);
        return dataset;
    }

    template<class Derived>
    template<class SpaceType>
    void ContinuousVector<Derived>::read(const H5DataSet<1>& dataset, const DataSpaceBase<SpaceType>& file_space) {
        const size_t length = dataset.getSize(0);
        resize(length);
        const auto mem_space = H5DataSpace<1>({length});
        dataset.read(data(), ScalarType::getH5DataType(), mem_space, file_space.asH5Type());
    }

    template<class Derived>
    H5DataSet<1> ContinuousVector<Derived>::write(H5Location& loc, const char* name) const {
        const auto space = H5DataSpace<1>({Base::getLength()});
        H5DataSet<1> dataset;
        if (loc.exists(name))
            dataset = loc.openDataSet<1>(name);
        else
            dataset = loc.createDataSet<1>(name, ScalarType::getH5DataType(), space);
        dataset.write(data(), ScalarType::getH5DataType(), space, space);
        return std::cref(dataset);
    }

    template<class Derived>
    template<class SpaceType>
    void ContinuousVector<Derived>::write(H5DataSet<1>& dataset, const DataSpaceBase<SpaceType>& file_space) const {
        const auto mem_space = H5DataSpace<1>({Base::getLength()});
        dataset.write(data(), ScalarType::getH5DataType(), mem_space, file_space.asH5Type());
    }
#endif
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
    inline void operator+=(ContinuousVector<Derived>& v1, const RValueVector<OtherDerived>& v2) {
        constexpr size_t Size1 = Traits<Derived>::SizeAtCompile;
        constexpr size_t Size2 = Traits<OtherDerived>::SizeAtCompile;
        static_assert(Size1 == Dynamic || Size2 == Dynamic || Size1 == Size2, "[Error]: Size mismatch between two vector");
        assert(v1.getLength() == v2.getLength());
        Internal::AddAssignImpl<Derived, OtherDerived, Internal::EnableSIMD<Derived, OtherDerived>::value>::run(v1, v2);
    }
}
