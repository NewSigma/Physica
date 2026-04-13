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
#include "../CompactVector.h"

namespace Physica {
    template<class Derived>
    auto CompactVector<Derived>::operator=(Scalar auto x) noexcept -> Derived& {
        if constexpr (SizeAtCompile == Dynamic) {
            if (x.isZero())
                zeros();
        }
        return Base::operator=(x);
    }

    template<class Derived>
    template<ExecutePolicy P>
    void CompactVector<Derived>::assign(Vector auto&& v) const noexcept {
        using V = std::remove_cvref<decltype(v)>::type;
        constexpr bool SameScalar = std::same_as<T, typename V::ScalarType>;
        constexpr bool Copyable = std::is_trivially_copyable<T>::value;
        if constexpr (V::isCompact() && SameScalar && Copyable) {
            if constexpr (isDiffable()) {
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
    template<int Size>
    auto CompactVector<Derived>::packet(size_t index) const noexcept -> SIMD<T, Size> {
        assert(index + Size <= Base::getLength());
        SIMD<T, Size> packet{};
        packet.load(Base::data_ptr(index));
        return packet;
    }

    template<class Derived>
    template<int Size>
    auto CompactVector<Derived>::packet(size_t index, size_t count) const noexcept -> SIMD<T, Size> {
        assert(index + count <= Base::getLength());
        assert(0 < count && count < Size && "[Error]: Invalid size for partial operation");
        SIMD<T, Size> packet{};
        packet.load(Base::data_ptr(index), count);
        return packet;
    }

    template<class Derived>
    void CompactVector<Derived>::writePacket(Packet auto packet, size_t index) noexcept {
        constexpr bool isSameScalar = std::is_same_v<ScalarType, typename Traits<decltype(packet)>::ScalarType>;
        if constexpr (isSameScalar)
            packet.store(Base::data_ptr(index));
        else
            Base::writePacket(packet, index);
    }

    template<class Derived>
    void CompactVector<Derived>::writePacket(Packet auto packet, size_t index, size_t count) noexcept {
        assert(index + count <= Base::getLength());
        assert(0 < count && count < packet.size() && "[Error]: Invalid size for partial operation");
        if constexpr (std::same_as<ScalarType, typename Traits<decltype(packet)>::ScalarType>)
            packet.store(Base::data_ptr(index), count);
        else
            Base::writePacket(packet, index, count);
    }

    template<class Derived>
    constexpr auto CompactVector<Derived>::view(this auto&& self) noexcept {
        return View<std::remove_reference_t<decltype(self)>>(self);
    }

    template<class Derived>
    template<size_t Length>
    auto CompactVector<Derived>::head(this auto&& self, size_t to) noexcept {
        using Self = decltype(self);
        return CompactVectorBlock<Self, Length>(std::forward<Self>(self), 0, to);
    }

    template<class Derived>
    template<size_t Length>
    auto CompactVector<Derived>::tail(this auto&& self, size_t from) noexcept {
        using Self = decltype(self);
        return CompactVectorBlock<Self, Length>(std::forward<Self>(self), from);
    }

    template<class Derived>
    template<size_t Length>
    auto CompactVector<Derived>::segment(this auto&& self, size_t from, size_t to) noexcept {
        using Self = decltype(self);
        return CompactVectorBlock<Self, Length>(std::forward<Self>(self), from, to);
    }

    template<class Derived>
    template<int Major, size_t Row, size_t Col>
    auto CompactVector<Derived>::reshape(this auto&& self, size_t row, size_t col) noexcept {
        using Self = decltype(self);
        return CompactReshapedVector<Self, Major, Row, Col>(std::forward<Self>(self), row, col);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto CompactVector<Derived>::reshape_row(this auto&& self, size_t row, size_t col) noexcept {
        return std::forward<decltype(self)>(self).template reshape<MatrixMajor::Row, Row, Col>(row, col);
    }

    template<class Derived>
    template<size_t Row, size_t Col>
    auto CompactVector<Derived>::reshape_col(this auto&& self, size_t row, size_t col) noexcept {
        return std::forward<decltype(self)>(self).template reshape<MatrixMajor::Col, Row, Col>(row, col);
    }

    template<class Derived>
    auto CompactVector<Derived>::reshape_like(this auto&& self, const Matrix auto& mat) noexcept {
        using M = std::remove_cvref_t<decltype(mat)>;
        constexpr auto Major = MatrixMajor::getMajor<M>();
        static_assert(Major != MatrixMajor::BothMajor, "[Error]: Cannot infer major from this matrix");
        return std::forward<decltype(self)>(self).template reshape<Major, M::RowAtCompile, M::ColAtCompile>(mat.getRow(), mat.getCol());
    }

    template<class Derived>
    auto CompactVector<Derived>::norm1() const noexcept -> CoDiff<Tr> {
        constexpr bool SmallVector = 0 < SizeAtCompile && SizeAtCompile <= 128;
        if constexpr (Internal::EnableMKL<Derived>::value && !SmallVector) {
            bool isSmallVector = Base::getLength() <= 128;
            return isSmallVector ? norm1_base() : norm1_mkl();
        }
        else
            return norm1_base();
    }

    template<class Derived>
    auto CompactVector<Derived>::norm1_base() const noexcept -> CoDiff<Tr> {
        return Base::norm1();
    }

    template<class Derived>
    auto CompactVector<Derived>::norm2() const noexcept -> CoDiff<Tr> {
        return norm2_base();
    }

    template<class Derived>
    auto CompactVector<Derived>::norm2_base() const noexcept -> CoDiff<Tr> {
        return Base::norm2();
    }
    /**
     * Prefer zeros() over simply assigning zeros for better performance.
     */
    template<class Derived>
    void CompactVector<Derived>::zeros() {
        if constexpr (Diffable<T>) {
            Base::getDerived().values().zeros();
            Base::getDerived().zero_grad();
        }
        else if constexpr (T::Prec != FloatMP)
            std::memset(data(), 0, Base::getLength() * sizeof(T));
        else
            Base::template operator=<T>(T(0));
    }
    /**
     * Read any Compact object and fetch enough scalars to fill self
     *
     * E.g. In optimization problems, we convert between complex vectors and real vectors.
     */
    template<class Derived>
    void CompactVector<Derived>::read(const auto& obj) noexcept {
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
    void CompactVector<Derived>::random_uniform() {
        if constexpr (R::MKL_Ready) {
            [[maybe_unused]] auto& gen = R::getInstance();
            [[maybe_unused]] const size_t length = Base::getLength() * (Base::isComplex() ? 2 : 1) * (Base::isForwardDiff() ? 2 : 1);
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
    void CompactVector<Derived>::random_normal() {
        if constexpr (R::MKL_Ready && !isForwardDiff()) {
            [[maybe_unused]] auto& gen = R::getInstance();
            [[maybe_unused]] const size_t length = Base::getLength() * (Base::isComplex() ? 2 : 1);
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
    auto CompactVector<Derived>::read(const H5Loc& loc, const char* name) -> const DataSetType {
        const auto dataset = loc.openDataSet<DataDim>(name);
        const size_t length = dataset.getSize(0);
        Base::resize(length);

        const auto memSpace = H5DataSpace<1>(length);
        if constexpr (isDiffable()) {
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
    auto CompactVector<Derived>::write(H5Loc& loc, const char* name) const -> DataSetType {
        const size_t length = Base::getLength();
        const auto memSpace = H5DataSpace<1>(length);
        DataSpaceType fileSpace;
        if constexpr (isDiffable())
            fileSpace = DataSpaceType({length, DiffOrder + 1});
        else
            fileSpace = memSpace;

        DataSetType dataset;
        if (loc.exists(name))
            dataset = loc.openDataSet<DataDim>(name);
        else
            dataset = loc.createDataSet<DataDim>(name, Tv::dtype_hdf5(), fileSpace);

        if constexpr (isDiffable()) {
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
    auto CompactVector<Derived>::data() noexcept {
        return Base::getDerived().data();
    }

    template<class Derived>
    auto CompactVector<Derived>::data() const noexcept {
        return Base::getDerived().data();
    }

    template<class Derived>
    auto CompactVector<Derived>::data_ptr(this auto&& self, size_t index) noexcept {
        return self.data() + index;
    }

    template<class Derived>
    __host__ __device__ consteval bool CompactVector<Derived>::isFastPacket() noexcept {
        return true;
    }
}
