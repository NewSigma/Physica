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

#include "Physica/Core/Exception/NoImplException.h"
#include "Physica/Core/Exception/MKL/VSL.h"

namespace Physica::Core {
    namespace Internal {
        template<LVector T1, Vector T2, bool enableSIMD>
        struct AddAssignImpl {
            static void run(T1& v1, const T2& v2) {
                using ScalarType = T1::ScalarType;
                for(size_t i = 0; i < v1.getLength(); ++i)
                    v1[i] += ScalarType(v2.calc(i));
            }
        };

        template<LVector T1, Vector T2>
        class AddAssignImpl<T1, T2, true> {
            constexpr static size_t Size1 = T1::SizeAtCompile;
            constexpr static size_t Size2 = T2::SizeAtCompile;
            constexpr static size_t SizeAtCompile = Size1 > Size2 ? Size1 : Size2;
            using PacketType = BestPacket<typename T1::ScalarType, SizeAtCompile>::Type;
        public:
            static void run(T1& v1, const T2& v2) {
                if constexpr (SizeAtCompile != Dynamic) {
                    constexpr size_t to = SizeAtCompile / PacketType::size() * PacketType::size();
                    for (size_t i = 0; i < to; i += PacketType::size()) {
                        const PacketType sum = v1.template packet<PacketType>(i) + v2.template packet<PacketType>(i);
                        v1.writePacket(i, sum);
                    }

                    constexpr size_t i = SizeAtCompile - SizeAtCompile % PacketType::size();
                    if constexpr (i != SizeAtCompile) {
                        constexpr size_t count = SizeAtCompile - i;
                        const PacketType sum = v1.template packetPartial<PacketType>(i, count) + v2.template packetPartial<PacketType>(i, count);
                        v1.writePacketPartial(i, count, sum);
                    }
                }
                else {
                    const size_t length = v1.getLength();
                    size_t i = 0;
                    const size_t to = length / PacketType::size() * PacketType::size();
                    for (; i < to; i += PacketType::size()) {
                        const PacketType sum = v1.template packet<PacketType>(i) + v2.template packet<PacketType>(i);
                        v1.writePacket(i, sum);
                    }
                    if (to != length) {
                        const size_t count = length - i;
                        const PacketType sum = v1.template packetPartial<PacketType>(i, count) + v2.template packetPartial<PacketType>(i, count);
                        v1.writePacketPartial(i, count, sum);
                    }
                }
            }
        };
    }

    template<class Derived>
    inline ContinuousVector<Derived>& ContinuousVector<Derived>::operator=(const This& v) {
        Base::operator=(v);
        return *this;
    }

    template<class Derived>
    inline ContinuousVector<Derived>& ContinuousVector<Derived>::operator=(This&& v) {
        return *this = v;
    }

    template<class Derived>
    template<Scalar T>
    inline Derived& ContinuousVector<Derived>::operator=(const T& x) {
        const size_t length = Base::getLength();
        if constexpr (isReverseDiff) {
            using TracerType = ScalarType::TracerType;
            TracerType::getInstance().reserve(length);
            for (size_t i = 0; i < length; ++i) {
                if constexpr (std::is_same<ScalarType, T>::value)
                    (*this)[i] = x.copy();
                else
                    (*this)[i] = ScalarType(x);
            }
        }
        else {
            const auto x1 = ScalarType(x);
            for (size_t i = 0; i < length; ++i)
                (*this)[i] = x1;
        }
        return Base::getDerived();
    }

    template<class Derived>
    template<Vector V>
    inline void ContinuousVector<Derived>::operator+=(const V& v) {
        constexpr size_t Size1 = SizeAtCompile;
        constexpr size_t Size2 = Traits<V>::SizeAtCompile;
        static_assert(Size1 == Dynamic || Size2 == Dynamic || Size1 == Size2, "[Error]: Size mismatch between two vector");
        assert(Base::getLength() == v.getLength());
        Internal::AddAssignImpl<Derived, V, Internal::EnableSIMD<Derived, V>::value>::run(Base::getDerived(), v);
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
            packet.load_partial(Base::data_ptr(index), count);
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
            Base::template writePacket<AnyPacket>(index, packet);
    }

    template<class Derived>
    template<class AnyPacket>
    inline void ContinuousVector<Derived>::writePacketPartial(size_t index, size_t count, const AnyPacket packet) {
        constexpr bool isSameScalar = std::is_same_v<ScalarType, typename Traits<AnyPacket>::ScalarType>;
        if constexpr (isSameScalar)
            packet.store_partial(Base::data_ptr(index), count);
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
    template<RandomGenerator R>
    inline void ContinuousVector<Derived>::random_uniform() {
        if constexpr (isReverseDiff) {
            using TracerType = ScalarType::TracerType;
            const size_t length = Base::getLength();
            TracerType::getInstance().reserve(length);
            for (size_t i = 0; i < length; ++i)
                this->operator[](i) = ScalarType::template random_uniform<R>();
        }
        else if constexpr (HasMKL()) {
            [[maybe_unused]] const size_t length = Base::getLength() * (Base::isComplex ? 2 : 1) * (Base::isForwardDiff ? 2 : 1);
            [[maybe_unused]] auto& gen = R::getInstance();
            if constexpr (ScalarType::Option == Float32)
                vslCheck(vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, gen, length, (float*)data(), 0, 1));
            else if constexpr (ScalarType::Option == Float64)
                vslCheck(vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, gen, length, (double*)data(), 0, 1));
            else
                Base::template random_uniform<R>();
        }
        else
            Base::template random_uniform<R>();
    }

    template<class Derived>
    template<RandomGenerator R>
    inline void ContinuousVector<Derived>::random_normal() {
        if constexpr (isReverseDiff) {
            using TracerType = ScalarType::TracerType;
            const size_t length = Base::getLength();
            TracerType::getInstance().reserve(length);
            for (size_t i = 0; i < length; ++i)
                this->operator[](i) = ScalarType::template random_normal<R>();
        }
        else if constexpr (HasMKL() && !isForwardDiff) {
            [[maybe_unused]] const size_t length = Base::getLength() * (Base::isComplex ? 2 : 1);
            [[maybe_unused]] auto& gen = R::getInstance();
            if constexpr (ScalarType::Option == Float32)
                vslCheck(vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, gen, length, (float*)data(), 0, 1));
            else if constexpr (ScalarType::Option == Float64)
                vslCheck(vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, gen, length, (double*)data(), 0, 1));
            else
                Base::template random_normal<R>();
        }
        else
            Base::template random_normal<R>();
    }

    template<class Derived>
    template<class Distribution, RandomGenerator R>
    inline void ContinuousVector<Derived>::random_any(Distribution& dist) {
        if constexpr (isReverseDiff) {
            using TracerType = ScalarType::TracerType;
            const size_t length = Base::getLength();
            TracerType::getInstance().reserve(length);
            for (size_t i = 0; i < length; ++i)
                this->operator[](i) = ScalarType::template random_any<decltype(dist), R>(dist);
        }
        else
            Base::template random_any<decltype(dist), R>(dist);
    }
#ifdef PHYSICA_HDF5
    template<class Derived>
    const ContinuousVector<Derived>::DataSetType ContinuousVector<Derived>::read(const H5Location& loc, const char* name) {
        const auto dataset = loc.openDataSet<DataDim>(name);
        const size_t length = dataset.getSize(0);
        resize(length);

        const auto memSpace = H5DataSpace<1>(length);
        if constexpr (isForwardDiff) {
            auto fileSpace = DataSpaceType({length, DiffOrder + 1});
            for (size_t i = 0; i <= DiffOrder; ++i) {
                fileSpace.selectHyperslab(H5S_SELECT_SET, {length, 1}, {0, i});
                dataset.read(data()[i], ValueType::getH5DataType(), memSpace, fileSpace);
            }
        }
        else
            dataset.read(data(), ValueType::getH5DataType(), memSpace, memSpace);
        return dataset;
    }

    template<class Derived>
    ContinuousVector<Derived>::DataSetType ContinuousVector<Derived>::write(H5Location& loc, const char* name) const {
        const size_t length = Base::getLength();
        const auto memSpace = H5DataSpace<1>(length);
        DataSpaceType fileSpace;
        if constexpr (isForwardDiff)
            fileSpace = DataSpaceType({length, DiffOrder + 1});
        else
            fileSpace = memSpace;

        DataSetType dataset;
        if (loc.exists(name))
            dataset = loc.openDataSet<DataDim>(name);
        else
            dataset = loc.createDataSet<DataDim>(name, ValueType::getH5DataType(), fileSpace);

        if constexpr (isForwardDiff) {
            for (size_t i = 0; i <= DiffOrder; ++i) {
                fileSpace.selectHyperslab(H5S_SELECT_SET, {length, 1}, {0, i});
                dataset.write(data()[i], ValueType::getH5DataType(), memSpace, fileSpace);
            }
        }
        else
            dataset.write(data(), ValueType::getH5DataType(), memSpace, fileSpace);
        return std::cref(dataset);
    }
#endif
    template<Vector T>
    std::ostream& operator<<(std::ostream& os, const ContinuousVector<T>& v) {
        using ScalarType = T::ScalarType;
        os.write(reinterpret_cast<const char*>(v.data_ptr(0)), v.getLength() * sizeof(ScalarType));
        return os;
    }

    template<Vector T>
    std::istream& operator>>(std::istream& is, ContinuousVector<T>& v) {
        using ScalarType = T::ScalarType;
        is.read(reinterpret_cast<char*>(v.data_ptr(0)), v.getLength() * sizeof(ScalarType));
        return is;
    }
}
