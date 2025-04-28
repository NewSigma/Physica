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

#include "Physica/Core/Exception/MKL/VSL.h"
#include "../ContinuousVector.h"

namespace Physica {
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
    template<Packet Pack>
    inline Pack ContinuousVector<Derived>::packet(size_t index) const {
        assert(index + Pack::size() <= Base::getLength());
        if constexpr (std::is_same_v<ScalarType, typename Traits<Pack>::ScalarType>) {
            Pack packet{};
            packet.load(Base::data_ptr(index));
            return packet;
        }
        else
            return Base::template packet<Pack>(index);
    }

    template<class Derived>
    template<Packet Pack>
    inline Pack ContinuousVector<Derived>::packetPartial(size_t index, size_t count) const {
        assert(index + count <= Base::getLength());
        if constexpr (std::is_same_v<ScalarType, typename Traits<Pack>::ScalarType>) {
            Pack packet{};
            packet.load_partial(Base::data_ptr(index), count);
            return packet;
        }
        else
            return Base::template packetPartial<Pack>(index, count);
    }

    template<class Derived>
    template<Packet Pack>
    inline void ContinuousVector<Derived>::writePacket(size_t index, const Pack packet) {
        constexpr bool isSameScalar = std::is_same_v<ScalarType, typename Traits<Pack>::ScalarType>;
        if constexpr (isSameScalar)
            packet.store(Base::data_ptr(index));
        else
            Base::template writePacket<Pack>(index, packet);
    }

    template<class Derived>
    template<Packet Pack>
    inline void ContinuousVector<Derived>::writePacketPartial(size_t index, size_t count, const Pack packet) {
        constexpr bool isSameScalar = std::is_same_v<ScalarType, typename Traits<Pack>::ScalarType>;
        if constexpr (isSameScalar)
            packet.store_partial(Base::data_ptr(index), count);
        else
            Base::template writePacketPartial<Pack>(index, count, packet);
    }

    template<class Derived>
    template<class U>
    void ContinuousVector<Derived>::reverse(const U& grad) const noexcept requires(isReverseDiff) {
        static_assert(std::same_as<typename ScalarType::GradType, typename U::ScalarType>, "[Error]: Inconsistent ScalarType");
        if constexpr (Scalar<U> || Vector<U>) {
            if constexpr (Vector<U>)
                assert(Base::getLength() == grad.getLength());
            Base::getConstCastDerived().grads() += grad;
        }
        else {
            static_assert(Matrix<U>, "[Error]: Unexpected type");
            assert(Base::getLength() == grad.getRow());
            reverse(grad.sum_cols());
        }
    }

    template<class Derived>
    template<size_t Length>
    inline auto ContinuousVector<Derived>::head(size_t to) noexcept {
        return BlockType<Length>(Base::getDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    inline const auto ContinuousVector<Derived>::head(size_t to) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    inline auto ContinuousVector<Derived>::tail(size_t from) noexcept {
        return BlockType<Length>(Base::getDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    inline const auto ContinuousVector<Derived>::tail(size_t from) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    inline auto ContinuousVector<Derived>::segment(size_t from, size_t to) noexcept {
        return BlockType<Length>(Base::getDerived(), from, to);
    }

    template<class Derived>
    template<size_t Length>
    inline const auto ContinuousVector<Derived>::segment(size_t from, size_t to) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from, to);
    }

    template<class Derived>
    void ContinuousVector<Derived>::zeros() {
        if constexpr (Diffable<T>) {
            Base::getDerived().values().zeros();
            Base::getDerived().grads().zeros();
        }
        else
            std::memset(data(), 0, Base::getLength() * sizeof(T));
    }

    template<class Derived>
    template<RNG R>
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
            if constexpr (ScalarType::Prec == Float32)
                check_vsl(vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, gen, length, (float*)data(), 0, 1));
            else if constexpr (ScalarType::Prec == Float64)
                check_vsl(vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, gen, length, (double*)data(), 0, 1));
            else
                Base::template random_uniform<R>();
        }
        else
            Base::template random_uniform<R>();
    }

    template<class Derived>
    template<RNG R>
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
            if constexpr (ScalarType::Prec == Float32)
                check_vsl(vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, gen, length, (float*)data(), 0, 1));
            else if constexpr (ScalarType::Prec == Float64)
                check_vsl(vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, gen, length, (double*)data(), 0, 1));
            else
                Base::template random_normal<R>();
        }
        else
            Base::template random_normal<R>();
    }

    template<class Derived>
    template<RNG R, class Distribution>
    inline void ContinuousVector<Derived>::random_any(Distribution& dist) {
        if constexpr (isReverseDiff) {
            using TracerType = ScalarType::TracerType;
            const size_t length = Base::getLength();
            TracerType::getInstance().reserve(length);
            for (size_t i = 0; i < length; ++i)
                this->operator[](i) = ScalarType::template random_any<R, decltype(dist)>(dist);
        }
        else
            Base::template random_any<R, decltype(dist)>(dist);
    }

#ifdef PHYSICA_HDF5
    template<class Derived>
    auto ContinuousVector<Derived>::read(const H5Loc& loc, const char* name) -> const DataSetType {
        const auto dataset = loc.openDataSet<DataDim>(name);
        const size_t length = dataset.getSize(0);
        resize(length);

        const auto memSpace = H5DataSpace<1>(length);
        if constexpr (isDiffable) {
            auto fileSpace = DataSpaceType({length, DiffOrder + 1});
            for (size_t i = 0; i <= DiffOrder; ++i) {
                fileSpace.selectHyperslab(H5S_SELECT_SET, {length, 1}, {0, i});
                dataset.read(data()[i], Tv::getH5DataType(), memSpace, fileSpace);
            }
        }
        else
            dataset.read(data(), Tv::getH5DataType(), memSpace, memSpace);
        return dataset;
    }

    template<class Derived>
    auto ContinuousVector<Derived>::write(H5Loc& loc, const char* name) const -> DataSetType {
        const size_t length = Base::getLength();
        const auto memSpace = H5DataSpace<1>(length);
        DataSpaceType fileSpace;
        if constexpr (isDiffable)
            fileSpace = DataSpaceType({length, DiffOrder + 1});
        else
            fileSpace = memSpace;

        DataSetType dataset;
        if (loc.exists(name))
            dataset = loc.openDataSet<DataDim>(name);
        else
            dataset = loc.createDataSet<DataDim>(name, Tv::getH5DataType(), fileSpace);

        if constexpr (isDiffable) {
            for (size_t i = 0; i <= DiffOrder; ++i) {
                fileSpace.selectHyperslab(H5S_SELECT_SET, {length, 1}, {0, i});
                dataset.write(data()[i], Tv::getH5DataType(), memSpace, fileSpace);
            }
        }
        else
            dataset.write(data(), Tv::getH5DataType(), memSpace, fileSpace);
        return std::cref(dataset);
    }
#endif

    template<class Derived>
    auto ContinuousVector<Derived>::values() noexcept -> ValuesRtnTy {
        return Base::getDerived();
    }

    template<class Derived>
    auto ContinuousVector<Derived>::values() const noexcept -> const ValuesRtnTy {
        return const_cast<This&>(*this).values();
    }

    template<class Derived>
    template<int GradOrder>
    auto ContinuousVector<Derived>::grads() const noexcept {
        static_assert(Diffable<ScalarType>, "[Error]: Undiffable vector does not have grads");
        return GradVector<ContinuousVector<Derived>, GradOrder>(Base::getDerived());
    }
}
