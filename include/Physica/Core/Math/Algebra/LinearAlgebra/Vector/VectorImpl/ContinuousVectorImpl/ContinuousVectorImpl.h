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

#include "Physica/Core/Exception/MKL/VSL.h"
#include "../ContinuousVector.h"

namespace Physica {
    template<class Derived>
    auto ContinuousVector<Derived>::operator=(Scalar auto x) noexcept -> Derived& {
        if constexpr (SizeAtCompile == Dynamic) {
            if (x.isZero())
                zeros();
        }
        return Base::operator=(x);
    }

    template<class Derived>
    template<ExecutePolicy P>
    void ContinuousVector<Derived>::assign(Vector auto&& v) const noexcept {
        using V = std::remove_cvref<decltype(v)>::type;
        constexpr bool SameScalar = std::same_as<T, typename V::ScalarType>;
        constexpr bool Copyable = std::is_trivially_copyable<T>::value;
        if constexpr (V::IsContinuous && SameScalar && Copyable) {
            if constexpr (isDiffable) {
                auto& x = Base::getDerived();
                x.values().template assign<P>(v.values());
                x.grads().template assign<P>(v.grads());
                return;
            }
            else {
                constexpr size_t Length = std::max(SizeAtCompile, V::SizeAtCompile);
                if constexpr (Length != Dynamic)
                    memcpy(v.data(), data(), Length * sizeof(T)); // Static memcpy
                else {
                    constexpr size_t Critical = 1024 * sizeof(float64); // Based on benchmark
                    constexpr bool MaybeBenefit = Length == Dynamic || (Length * sizeof(T) > Critical);
                    size_t size = Base::getLength() * sizeof(T);
                    if (MaybeBenefit && (size > Critical))
                        memcpy(v.data(), data(), size); // Reuse size, because v.read(data()) shows IR regression
                    else
                        Base::template assign_base<P>(v);
                }
            }
        }
        else
            Base::template assign_base<P>(v);
    }

    template<class Derived>
    template<Packet Pack>
    Pack ContinuousVector<Derived>::packet(size_t index) const {
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
    Pack ContinuousVector<Derived>::packetPartial(size_t index, size_t count) const {
        assert(index + count <= Base::getLength());
        assert(0 < count && count < Pack::size() && "[Error]: Invalid size for partial operation");
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
    void ContinuousVector<Derived>::writePacket(size_t index, const Pack packet) {
        constexpr bool isSameScalar = std::is_same_v<ScalarType, typename Traits<Pack>::ScalarType>;
        if constexpr (isSameScalar)
            packet.store(Base::data_ptr(index));
        else
            Base::template writePacket<Pack>(index, packet);
    }

    template<class Derived>
    template<Packet Pack>
    void ContinuousVector<Derived>::writePacketPartial(size_t index, size_t count, const Pack packet) {
        assert(index + count <= Base::getLength());
        assert(0 < count && count < Pack::size() && "[Error]: Invalid size for partial operation");
        if constexpr (std::same_as<ScalarType, typename Traits<Pack>::ScalarType>)
            packet.store_partial(Base::data_ptr(index), count);
        else
            Base::template writePacketPartial<Pack>(index, count, packet);
    }

    template<class Derived>
    template<size_t Length>
    auto ContinuousVector<Derived>::head(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return ContinuousVectorBlock<Self, Length>(std::forward<Self>(self), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    auto ContinuousVector<Derived>::tail(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return ContinuousVectorBlock<Self, Length>(std::forward<Self>(self), from);
    }

    template<class Derived>
    template<size_t Length>
    auto ContinuousVector<Derived>::segment(this auto&& self, size_t from, size_t to) noexcept {
        using Self = decltype(self);
        return ContinuousVectorBlock<Self, Length>(std::forward<Self>(self), from, to);
    }

    template<class Derived>
    auto ContinuousVector<Derived>::norm1() const noexcept -> CoDiff<Tr> {
        constexpr bool SmallVector = 0 < SizeAtCompile && SizeAtCompile <= 128;
        if constexpr (Internal::EnableMKL<Derived>::value && !SmallVector) {
            bool isSmallVector = Base::getLength() <= 128;
            return isSmallVector ? norm1_base() : norm1_mkl();
        }
        else
            return norm1_base();
    }

    template<class Derived>
    auto ContinuousVector<Derived>::norm1_base() const noexcept -> CoDiff<Tr> {
        return Base::norm1();
    }

    template<class Derived>
    auto ContinuousVector<Derived>::norm2() const noexcept -> CoDiff<Tr> {
        return norm2_base();
    }

    template<class Derived>
    auto ContinuousVector<Derived>::norm2_base() const noexcept -> CoDiff<Tr> {
        return Base::norm2();
    }
    /**
     * Prefer zeros() over simply assigning zeros for better performance.
     */
    template<class Derived>
    void ContinuousVector<Derived>::zeros() {
        if constexpr (Diffable<T>) {
            Base::getDerived().values().zeros();
            Base::getDerived().grads().zeros();
        }
        else if constexpr (T::Prec != FloatMP)
            std::memset(data(), 0, Base::getLength() * sizeof(T));
        else
            Base::template operator=<T>(T(0));
    }
    /**
     * Read any continuous object and fetch enough scalars to fill self
     *
     * E.g. In optimization problems, we convert between complex vectors and real vectors.
     */
    template<class Derived>
    void ContinuousVector<Derived>::read(const auto& obj) noexcept {
        using O = decltype(obj);
        if constexpr (Vector<O>) {
            using U = std::remove_cvref_t<O>::ScalarType;
            static_assert(T::Prec == U::Prec);
            size_t size = Base::getLength() * sizeof(T);
            assert(size <= obj.getLength() * sizeof(U));
            memcpy(data(), obj.data(), size);
        }
        else {
            static_assert(Matrix<O>, "[Error]: Unexpected type");
            read(obj.flatten());
        }
    }

    template<class Derived>
    template<RNG R>
    void ContinuousVector<Derived>::random_uniform() {
        if constexpr (R::MKL_Ready) {
            [[maybe_unused]] auto& gen = R::getInstance();
            [[maybe_unused]] const size_t length = Base::getLength() * (Base::isComplex ? 2 : 1) * (Base::isForwardDiff ? 2 : 1);
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
    void ContinuousVector<Derived>::random_normal() {
        if constexpr (R::MKL_Ready && !isForwardDiff) {
            [[maybe_unused]] auto& gen = R::getInstance();
            [[maybe_unused]] const size_t length = Base::getLength() * (Base::isComplex ? 2 : 1);
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

#ifdef PHYSICA_HDF5
    template<class Derived>
    auto ContinuousVector<Derived>::read(const H5Loc& loc, const char* name) -> const DataSetType {
        const auto dataset = loc.openDataSet<DataDim>(name);
        const size_t length = dataset.getSize(0);
        Base::resize(length);

        const auto memSpace = H5DataSpace<1>(length);
        if constexpr (isDiffable) {
            auto fileSpace = DataSpaceType({length, DiffOrder + 1});
            for (size_t i = 0; i <= DiffOrder; ++i) {
                fileSpace.selectHyperslab(H5S_SELECT_SET, {length, 1}, {0, i});
                dataset.read(data().get(i), Tv::dtype_hdf5(), memSpace, fileSpace);
            }
        }
        else
            dataset.read(data(), Tv::dtype_hdf5(), memSpace, memSpace);
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
            dataset = loc.createDataSet<DataDim>(name, Tv::dtype_hdf5(), fileSpace);

        if constexpr (isDiffable) {
            for (size_t i = 0; i <= DiffOrder; ++i) {
                fileSpace.selectHyperslab(H5S_SELECT_SET, {length, 1}, {0, i});
                dataset.write(data().get(i), Tv::dtype_hdf5(), memSpace, fileSpace);
            }
        }
        else
            dataset.write(data(), Tv::dtype_hdf5(), memSpace, fileSpace);
        return std::cref(dataset);
    }
#endif

    template<class Derived>
    auto ContinuousVector<Derived>::data() noexcept {
        return Base::getDerived().data();
    }

    template<class Derived>
    auto ContinuousVector<Derived>::data() const noexcept {
        return Base::getDerived().data();
    }

    template<class Derived>
    auto ContinuousVector<Derived>::data_ptr(this auto&& self, size_t index) noexcept {
        return self.data() + index;
    }
}
