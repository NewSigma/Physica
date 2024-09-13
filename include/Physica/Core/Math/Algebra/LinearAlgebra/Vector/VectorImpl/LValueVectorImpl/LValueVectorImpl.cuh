/*
 * Copyright 2023-2024 Weibo He.
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

namespace Physica::Core {
    namespace Internal {
        template<class Derived, class ScalarType>
        __global__ void assignConst_kernel(
                Physica::PlainStruct<Derived> target,
                ScalarType constant) {
            using namespace Physica::Core;
            const unsigned int delta = gridDim.x * blockDim.x;
            const unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
            const size_t length = target.getDerived().getLength();
            for (unsigned int shift = 0; shift < length; shift += delta) {
                const unsigned int index = id + shift;
                if (index < length)
                    target.getDerived()[index] = constant;
            }
        }
    }

    template<class Derived>
    __host__ __device__
    inline device_obj<LValueVector<Derived>>&
    device_obj<LValueVector<Derived>>::operator=(const device_obj<LValueVector<Derived>>& obj) {
        Base::operator=(obj);
        return *this;
    }

    template<class Derived>
    __host__ __device__
    inline device_obj<LValueVector<Derived>>&
    device_obj<LValueVector<Derived>>::operator=(device_obj<LValueVector<Derived>>&& obj) noexcept {
        Base::operator=(std::move(obj));
        return *this;
    }

    template<class Derived>
    template<class OtherDerived>
    __host__ __device__
    device_obj<Derived>& device_obj<LValueVector<Derived>>::operator=(const device_obj<RValueVector<OtherDerived>>& v) {
        if constexpr (IsHost())
            Base::getDerived().resize(v.getLength());
        v.assignTo(*this);
        return Base::getDerived();
    }

    template<class Derived>
    template<class AnyScalar>
    inline device_obj<Derived>& device_obj<LValueVector<Derived>>::operator=(const ScalarBase<AnyScalar>& s) {
        using namespace Physica;
        constexpr unsigned int WarpSize = Physica::CUDADevAttr::WarpSize;
        const int numBlock = (Base::getLength() + WarpSize - 1) / WarpSize;
        const int numThread = WarpSize;
        Internal::assignConst_kernel<<<numBlock, numThread, 0, CUDAContext::getInstance()>>>(asStruct(Base::getDerived()), s.getDerived());
        check(cudaGetLastError());
        return Base::getDerived();
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline device_obj<LVectorBlock<Derived, Length>> device_obj<LValueVector<Derived>>::head(size_t to) {
        return {Base::getDerived(), 0, to};
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline const device_obj<LVectorBlock<Derived, Length>> device_obj<LValueVector<Derived>>::head(size_t to) const {
        return {Base::getConstCastDerived(), 0, to};
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline device_obj<LVectorBlock<Derived, Length>> device_obj<LValueVector<Derived>>::tail(size_t from) {
        return {Base::getDerived(), from};
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline const device_obj<LVectorBlock<Derived, Length>> device_obj<LValueVector<Derived>>::tail(size_t from) const {
        return {Base::getConstCastDerived(), from};
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline device_obj<LVectorBlock<Derived, Length>> device_obj<LValueVector<Derived>>::segment(size_t from, size_t to) {
        return {Base::getDerived(), from, to};
    }

    template<class Derived>
    template<size_t Length>
    __host__ __device__ inline const device_obj<LVectorBlock<Derived, Length>> device_obj<LValueVector<Derived>>::segment(size_t from, size_t to) const {
        return {Base::getConstCastDerived(), from, to};
    }

    template<class Derived>
    __host__ __device__ inline typename device_obj<LValueVector<Derived>>::ScalarType*
    device_obj<LValueVector<Derived>>::data_ptr(size_t index) {
        return Base::getDerived().data_ptr(index);
    }
    
    template<class Derived>
    __host__ __device__ inline const typename device_obj<LValueVector<Derived>>::ScalarType*
    device_obj<LValueVector<Derived>>::data_ptr(size_t index) const {
        return Base::getDerived().data_ptr(index);
    }
}