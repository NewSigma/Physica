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

#include "ContinuousVectorImpl/ContinuousVectorBlock.h"

namespace Physica {
    template<class Derived>
    class ContinuousVector : public LValueVector<Derived> {
        using Base = LValueVector<Derived>;
        using This = ContinuousVector<Derived>;
        template<size_t Length>
        using BlockType = ContinuousVectorBlock<Derived, Length>;
    public:
        using typename Base::ScalarType;
        using typename Base::PacketType;
        using Base::SizeAtCompile;
        using Base::isForwardDiff;
        using Base::isReverseDiff;
        using Base::isDiffable;
    protected:
        using typename Base::T;
        using typename Base::Tr;
        using typename Base::Tv;
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    private:
        constexpr static int DiffOrder = ScalarType::Order;
        constexpr static int DataDim = 1 + (DiffOrder > 0);
        using DataSetType = H5DataSet<DataDim>;
        using DataSpaceType = H5DataSpace<DataDim>;

        using IterF = PtrIteratorF<Derived>;
        using IterR = PtrIteratorR<Derived>;
        using IterCF = PtrIteratorF<const Derived>;
        using IterCR = PtrIteratorR<const Derived>;

        using ValuesRtnTy = std::conditional<Diffable<ScalarType>, ValueVector<This>, Derived&>::type;
    public:
        ~ContinuousVector() = default;
        /* Operators */
        This& operator=(const This& v) = delete;
        This& operator=(This&& v) noexcept = delete;
        Derived& operator=(const Scalar auto& x);
        template<ExecutePolicy P = Sequential>
        Derived& operator=(const Vector auto& v);
        using Base::operator+=;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto&& v) const noexcept;
        template<ExecutePolicy P = Sequential>
        void assign_base(Vector auto& v) const noexcept;
        void assign_mkl(Vector auto& v) const noexcept;

        template<Packet Pack> [[nodiscard]] Pack packet(size_t index) const;
        template<Packet Pack> [[nodiscard]] Pack packetPartial(size_t index, size_t count) const;
        template<Packet Pack> void writePacket(size_t index, Pack packet);
        template<Packet Pack> void writePacketPartial(size_t index, size_t count, Pack packet);

        void reverse(const auto& grad) const noexcept requires(isReverseDiff);

        void resize(const Vector auto& x) { resize(x.getLength()); }
        auto resize(size_t length) { return Base::getDerived().resize(length); }
        template<Vector T> void toDevice(device_obj<ContinuousVector<T>>& obj) const;
        template<Vector T> void toDeviceAsync(device_obj<ContinuousVector<T>>& obj) const;
        [[nodiscard]] auto toNumpy() const;

        template<size_t Length = Dynamic>
        [[nodiscard]] auto head(size_t to = Length) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] const auto head(size_t to = Length) const noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] auto tail(size_t from) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] const auto tail(size_t from) const noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] auto segment(size_t from, size_t to) noexcept;
        template<size_t Length = Dynamic>
        [[nodiscard]] const auto segment(size_t from, size_t to) const noexcept;

        [[nodiscard]] auto begin() noexcept { return IterF(data()); }
        [[nodiscard]] auto begin() const noexcept { return cbegin(); }
        [[nodiscard]] auto cbegin() const noexcept { return IterCF(data()); }
        [[nodiscard]] auto end() noexcept { return IterF(data() + Base::getLength()); }
        [[nodiscard]] auto end() const noexcept { return cend(); }
        [[nodiscard]] auto cend() const noexcept { return IterCF(data() + Base::getLength()); }
        [[nodiscard]] auto rbegin() noexcept { return IterR(data() + Base::getLength() - 1); }
        [[nodiscard]] auto rbegin() const noexcept { return crbegin(); }
        [[nodiscard]] auto crbegin() const noexcept { return IterCR(data() + Base::getLength() - 1); }
        [[nodiscard]] auto rend() noexcept { return IterR(data() - 1); }
        [[nodiscard]] auto rend() const noexcept { return crend(); }
        [[nodiscard]] auto crend() const noexcept { return IterCR(data() - 1); }

        [[nodiscard]] CoDiff<Tr> norm1() const noexcept;
        [[nodiscard]] CoDiff<Tr> norm1_base() const noexcept;
        [[nodiscard]] Tr norm1_mkl() const noexcept;
        [[nodiscard]] CoDiff<Tr> norm2() const noexcept;
        [[nodiscard]] CoDiff<Tr> norm2_base() const noexcept;
        [[nodiscard]] Tr norm2_mkl() const noexcept;

        void zeros();
        template<RNG R>
        void random_uniform();
        template<RNG R>
        void random_normal();
        template<RNG R>
        void random_any(auto& distribution);

        const DataSetType read(const H5Loc& loc, const char* name);
        DataSetType write(H5Loc& loc, const char* name) const;

        [[nodiscard]] ValuesRtnTy values() noexcept;
        [[nodiscard]] const ValuesRtnTy values() const noexcept;
        template<int GradOrder = 1>
        auto grads() const noexcept;
        /* Getters */
        [[nodiscard]] PtrTy data() { return Base::data_ptr(0); }
        [[nodiscard]] ConstPtrTy data() const { return Base::data_ptr(0); }
    protected:
        ContinuousVector() = default;
        ContinuousVector(const This&) = default;
        ContinuousVector(This&&) noexcept = default;
        /* Friends */
        friend class device_obj<ContinuousVector<Derived>>;
    };
}

#include "ContinuousVectorImpl/ContinuousVectorImpl.h"
#ifdef PHYSICA_MKL
    #include "ContinuousVectorImpl/ContinuousVector_MKL.h"
#endif
#include "ContinuousVectorImpl/VectorConvert.h"
