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
    Derived& ContinuousVector<Derived>::operator=(const Scalar auto& x) {
        if constexpr (SizeAtCompile == Dynamic) {
            if (x.isZero())
                zeros();
        }
        return Base::operator=(x);
    }

    template<class Derived>
    template<ExecutePolicy P>
    Derived& ContinuousVector<Derived>::operator=(const Vector auto& v) {
        return Base::operator=(v);
    }

    template<class Derived>
    template<ExecutePolicy P>
    void ContinuousVector<Derived>::assign(Vector auto& v) const noexcept {
        using V = std::remove_cvref<decltype(v)>::type;
        constexpr size_t Length = std::max(SizeAtCompile, V::SizeAtCompile);
        constexpr size_t Critical = 1024 * sizeof(float64); // Based on benchmark
        constexpr bool Beneficial = Length == Dynamic || (Length * sizeof(T) > Critical);
        constexpr bool isContinuous = is_continuous<V>::value;
        constexpr bool SameScalar = std::same_as<T, typename V::ScalarType>;
        constexpr bool Copyable = std::is_trivially_copyable<T>::value;
        if constexpr (Beneficial && isContinuous && SameScalar && Copyable) {
            if (Base::getLength() * sizeof(T) > Critical) {
                if constexpr (isDiffable) {
                    Base::getDerived().values().template assign<P>(v.values());
                    Base::getDerived().grads().template assign<P>(v.grads());
                }
                else
                    memcpy(v.data(), data(), Base::getLength() * sizeof(T));
            }
            else
                assign_base<P>(v);
        }
        else
            assign_base<P>(v);
    }

    template<class Derived>
    template<ExecutePolicy P>
    void ContinuousVector<Derived>::assign_base(Vector auto& v) const noexcept {
        Base::template assign<P>(v);
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
        constexpr bool isSameScalar = std::is_same_v<ScalarType, typename Traits<Pack>::ScalarType>;
        if constexpr (isSameScalar)
            packet.store_partial(Base::data_ptr(index), count);
        else
            Base::template writePacketPartial<Pack>(index, count, packet);
    }

    template<class Derived>
    void ContinuousVector<Derived>::reverse(const auto& grad) const noexcept requires(isReverseDiff) {
        using U = std::remove_cvref_t<decltype(grad)>;
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
    auto ContinuousVector<Derived>::head(size_t to) noexcept {
        return BlockType<Length>(Base::getDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    const auto ContinuousVector<Derived>::head(size_t to) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    auto ContinuousVector<Derived>::tail(size_t from) noexcept {
        return BlockType<Length>(Base::getDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    const auto ContinuousVector<Derived>::tail(size_t from) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from);
    }

    template<class Derived>
    template<size_t Length>
    auto ContinuousVector<Derived>::segment(size_t from, size_t to) noexcept {
        return BlockType<Length>(Base::getDerived(), from, to);
    }

    template<class Derived>
    template<size_t Length>
    const auto ContinuousVector<Derived>::segment(size_t from, size_t to) const noexcept {
        return BlockType<Length>(Base::getConstCastDerived(), from, to);
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
    auto ContinuousVector<Derived>::norm2() const noexcept -> CoDiff<Tr> {
        return norm2_base();
    }

    template<class Derived>
    auto ContinuousVector<Derived>::norm2_base() const noexcept -> CoDiff<Tr> {
        return Base::norm2();
    }

    template<class Derived>
    auto ContinuousVector<Derived>::norm1_base() const noexcept -> CoDiff<Tr> {
        return Base::norm1();
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

    template<class Derived>
    template<RNG R>
    void ContinuousVector<Derived>::random_uniform() {
        if constexpr (R::MKL_Ready) {
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
    void ContinuousVector<Derived>::random_normal() {
        if constexpr (R::MKL_Ready && !isForwardDiff) {
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
    template<RNG R>
    void ContinuousVector<Derived>::random_any(auto& distribution) {
        Base::template random_any<R>(distribution);
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
                dataset.read(data()[i], Tv::dtype_hdf5(), memSpace, fileSpace);
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
                dataset.write(data()[i], Tv::dtype_hdf5(), memSpace, fileSpace);
            }
        }
        else
            dataset.write(data(), Tv::dtype_hdf5(), memSpace, fileSpace);
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
