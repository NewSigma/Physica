/*
 * Copyright 2023 WeiBo He.
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

#include "ContinuousVector.h"
#include "LValueVector.cuh"
#include "ContinuousVectorBlock.cuh"

namespace Physica::Core {
    template<class Derived>
    class device_obj<ContinuousVector<Derived>> : public device_obj<LValueVector<Derived>> {
        using Base = device_obj<LValueVector<Derived>>;
        template<size_t Length> using BlockType = device_obj<ContinuousVectorBlock<Derived, Length>>;
    public:
        using host_obj = ContinuousVector<Derived>;
        using typename Base::ScalarType;
    public:
        ~device_obj() = default;
        /* Operators */
        inline device_obj& operator=(const device_obj& obj);
        inline device_obj& operator=(device_obj&& obj) noexcept;
        using Base::operator=;
        /* Operations */
        void resize(size_t length) { Base::getDerived().resize(length); }
        template<class OtherDerived>
        void toHost(ContinuousVector<OtherDerived>& obj) const;

        template<size_t Length = Dynamic> __host__ __device__ inline BlockType<Length> head(size_t to);
        template<size_t Length = Dynamic> __host__ __device__ inline const BlockType<Length> head(size_t to) const;
        template<size_t Length = Dynamic> __host__ __device__ inline BlockType<Length> tail(size_t from);
        template<size_t Length = Dynamic> __host__ __device__ inline const BlockType<Length> tail(size_t from) const;
        template<size_t Length = Dynamic> __host__ __device__ inline BlockType<Length> segment(size_t from, size_t to);
        template<size_t Length = Dynamic> __host__ __device__ inline const BlockType<Length> segment(size_t from, size_t to) const;
        /* Getters */
        [[nodiscard]] __device__ ScalarType* data() { return Base::getDerived().data(); }
        [[nodiscard]] __host__ __device__ const ScalarType* data() const { return Base::getDerived().data(); }
    protected:
        device_obj() = default;
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
    };

    template<class Derived>
    inline device_obj<ContinuousVector<Derived>>&
    device_obj<ContinuousVector<Derived>>::operator=(const device_obj<ContinuousVector<Derived>>& obj) {
        Base::operator=(obj);
        return *this;
    }
    
    template<class Derived>
    inline device_obj<ContinuousVector<Derived>>&
    device_obj<ContinuousVector<Derived>>::operator=(device_obj<ContinuousVector<Derived>>&& obj) noexcept {
        Base::operator=(std::move(obj));
        return *this;
    }

    template<class Derived>
    template<class OtherDerived>
    void device_obj<ContinuousVector<Derived>>::toHost(ContinuousVector<OtherDerived>& obj) const {
        obj.resize(Base::getLength());
        cudaCheck(cudaMemcpy(obj.data(), data(), Base::getLength() * sizeof(ScalarType), cudaMemcpyKind::cudaMemcpyDeviceToHost));
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline typename device_obj<ContinuousVector<Derived>>::template BlockType<Length>
    device_obj<ContinuousVector<Derived>>::head(size_t to) {
        return {Base::getDerived(), 0, to};
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline const typename device_obj<ContinuousVector<Derived>>::template BlockType<Length>
    device_obj<ContinuousVector<Derived>>::head(size_t to) const {
        return {Base::getConstCastDerived(), 0, to};
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline typename device_obj<ContinuousVector<Derived>>::template BlockType<Length>
    device_obj<ContinuousVector<Derived>>::tail(size_t from) {
        return {Base::getDerived(), from};
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline const typename device_obj<ContinuousVector<Derived>>::template BlockType<Length>
    device_obj<ContinuousVector<Derived>>::tail(size_t from) const {
        return {Base::getConstCastDerived(), from};
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline typename device_obj<ContinuousVector<Derived>>::template BlockType<Length>
    device_obj<ContinuousVector<Derived>>::segment(size_t from, size_t to) {
        return {Base::getDerived(), from, to};
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline const typename device_obj<ContinuousVector<Derived>>::template BlockType<Length>
    device_obj<ContinuousVector<Derived>>::segment(size_t from, size_t to) const {
        return {Base::getConstCastDerived(), from, to};
    }

    template<class Derived>
    template<class OtherDerived>
    void ContinuousVector<Derived>::toDevice(device_obj<ContinuousVector<OtherDerived>>& obj) const {
        obj.resize(Base::getLength());
        const device_obj<ContinuousVector<OtherDerived>>& const_obj = obj;
        cudaCheck(cudaMemcpy((void*)const_obj.data(), data(), Base::getLength() * sizeof(ScalarType), cudaMemcpyKind::cudaMemcpyHostToDevice));
    }
}
