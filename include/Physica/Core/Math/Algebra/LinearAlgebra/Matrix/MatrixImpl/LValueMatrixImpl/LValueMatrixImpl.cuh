/*
 * Copyright 2024-2025 Weibo He.
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

#include "../LValueMatrix.cuh"

namespace Physica {
    template<class Derived>
    template<Matrix M>
    __host__ __device__ device_obj<Derived>& device_obj<LValueMatrix<Derived>>::operator=(const M& m) requires(CUDA<M>) {
        static_assert(RowAtCompile == Dynamic || M::RowAtCompile == Dynamic || RowAtCompile == M::RowAtCompile, "[Error]: Incompatible row number");
        static_assert(ColAtCompile == Dynamic || M::ColAtCompile == Dynamic || ColAtCompile == M::ColAtCompile, "[Error]: Incompatible col number");
        auto& target = Base::getDerived();
        target.resize(m.getRow(), m.getCol());
        m.assign(target);
        return target;
    }

    template<class Derived>
    template<Scalar T>
    __host__ __device__ device_obj<Derived>& device_obj<LValueMatrix<Derived>>::operator=(const T& x) {
        if (IsHost()) {
            const auto config = Base::makeKernelConfig();
            CUDAExecutor::launch([m_ = asStruct(Base::getDerived()), x] __device__() mutable {
                const auto& m = m_.getDerived();
                const size_t maxMinor = m.getMaxMinor();
                const size_t major = blockIdx.y;
                const size_t minor = blockIdx.x * blockDim.x + threadIdx.x;
                if (minor < maxMinor)
                    m.refFromMajorMinor(major, minor) = x;
            }, config.first, config.second);
        }
        else {
            for (size_t i = 0; i < Base::getMaxMajor(); ++i)
                for (size_t j = 0; j < Base::getMaxMinor(); ++j)
                    refFromMajorMinor(i, j) = x;
            return Base::getDerived();
        }
    }

    template<class Derived>
    __device__ auto device_obj<LValueMatrix<Derived>>::operator()(size_t row, size_t col) -> RefTy {
        return *data_ptr(row, col);
    }

    template<class Derived>
    __device__ auto device_obj<LValueMatrix<Derived>>::operator()(size_t row, size_t col) const -> ConstRefTy {
        return *data_ptr(row, col);
    }

    template<class Derived>
    template<Matrix M>
    void device_obj<LValueMatrix<Derived>>::reverse(const M& grad) const noexcept requires(isReverseDiff) {
        static_assert(std::same_as<typename ScalarType::GradType, typename M::ScalarType>, "[Error]: Inconsistent ScalarType");
        assert(Base::getRow() == grad.getRow());
        assert(Base::getCol() == grad.getCol());
        Base::getConstCastDerived().grads() += grad;
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::flatten() {
        return device_obj<FlattenL<Derived>>(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ const auto device_obj<LValueMatrix<Derived>>::flatten() const {
        return device_obj<FlattenL<Derived>>(Base::getDerived());
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::data_ptr(size_t row, size_t col) -> PtrTy {
        assert(row < Base::getRow());
        assert(col < Base::getCol());
        return Base::getDerived().data_ptr(row, col);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<LValueMatrix<Derived>>::data_ptr(size_t row, size_t col) const -> ConstPtrTy {
        return const_cast<This&>(*this).data_ptr(row, col);
    }

    template<class Derived>
    __device__ inline auto device_obj<LValueMatrix<Derived>>::refFromMajorMinor(size_t major, size_t minor) -> RefTy {
        assert(major < Base::getDerived().getMaxMajor());
        assert(minor < Base::getDerived().getMaxMinor());
        const size_t r = MatrixOption::rowFromMajorMinor<Derived>(major, minor);
        const size_t c = MatrixOption::colFromMajorMinor<Derived>(major, minor);
        return Base::getDerived()(r, c);
    }

    template<class Derived>
    __device__ inline auto device_obj<LValueMatrix<Derived>>::refFromMajorMinor(size_t major, size_t minor) const -> ConstRefTy {
        return const_cast<This&>(*this).refFromMajorMinor(major, minor);
    }
}
