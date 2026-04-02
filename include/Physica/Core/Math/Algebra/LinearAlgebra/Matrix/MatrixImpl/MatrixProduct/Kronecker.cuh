/*
 * Copyright 2025-2026 Weibo He.
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

#include "../RValueMatrix.cuh"

namespace Physica {
    template<Matrix M1, Matrix M2>
    class device_obj<Kronecker<M1, M2>> : public device_obj<RValueMatrix<Kronecker<M1, M2>>> {
        using host_obj = Kronecker<M1, M2>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using Ref1 = add_device_obj<M1>::type;
        using Ref2 = add_device_obj<M2>::type;
    protected:
        using typename Base::T;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M1>>> m1;
        PlainStruct<add_device_obj_t<std::remove_reference_t<M2>>> m2;
    public:
        __host__ __device__ device_obj(Ref1 m1, Ref2 m2);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto&& target) const;
        void assign_add(Matrix auto&& target) const;

        [[nodiscard]] __device__ T calc(size_t row, size_t col) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept;
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept;
        [[nodiscard]] __host__ __device__ auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto&& getRHS(this auto&&) noexcept;
    private:
        void assign_identity(Matrix auto&& target) const;
        void assign_add_identity(Matrix auto&& target) const;
        void assign_diagonal(Matrix auto&& target) const;
        void assign_add_diagonal(Matrix auto&& target) const;
    };

    template<Matrix M1, Matrix M2>
    __host__ __device__ device_obj<Kronecker<M1, M2>>::device_obj(Ref1 m1, Ref2 m2) : m1(asStruct(m1)), m2(asStruct(m2)) {}

    template<Matrix M1, Matrix M2>
    void device_obj<Kronecker<M1, M2>>::assign(Matrix auto&& target) const {
        if constexpr (instanceof_tx<IdentityMatrix, M1>)
            assign_identity(target);
        else if constexpr (instanceof_tx<DiagMatrix, M1>)
            assign_diagonal(target);
        else {
            target.zeros();
            auto assign_general = [lhs_ = m1, rhs_ = m2, target_ = asStruct(target)] __device__() mutable {
                const auto& lhs = lhs_.getDerived();
                const auto& rhs = rhs_.getDerived();
                unsigned int r = blockIdx.y;
                unsigned int c = blockIdx.z;
                unsigned int offsetR = r * rhs.getRow();
                unsigned int offsetC = c * rhs.getCol();
                (lhs.calc(r, c) * rhs).assign(target_.getDerived().block(offsetR, rhs.getRow(), offsetC, rhs.getCol()), ThreadBlock{});
            };

            const auto& lhs = getLHS();
            const auto& rhs = getRHS();
            const auto numThread = std::min((unsigned int)rhs.getSize(), CUDADevAttr::DefaultThreadsPerBlock);
            const unsigned int numBlockX = 1;
            const unsigned int numBlockY = lhs.getRow();
            const unsigned int numBlockZ = lhs.getCol();
            CUDAExecutor::launch(assign_general, KernelConfig({numBlockX, numBlockY, numBlockZ}, numThread));
        }
    }

    template<Matrix M1, Matrix M2>
    void device_obj<Kronecker<M1, M2>>::assign_add(Matrix auto&& target) const {
        if constexpr (instanceof_tx<IdentityMatrix, M1>)
            assign_add_identity(target);
        else if constexpr (instanceof_tx<DiagMatrix, M1>)
            assign_add_diagonal(target);
        else {
            auto assign_add_general = [lhs_ = m1, rhs_ = m2, target_ = asStruct(target)] __device__() mutable {
                const auto& lhs = lhs_.getDerived();
                const auto& rhs = rhs_.getDerived();
                unsigned int r = blockIdx.y;
                unsigned int c = blockIdx.z;
                unsigned int offsetR = r * rhs.getRow();
                unsigned int offsetC = c * rhs.getCol();
                (lhs.calc(r, c) * rhs).assign_add(target_.getDerived().block(offsetR, rhs.getRow(), offsetC, rhs.getCol()), ThreadBlock{});
            };

            const auto& lhs = getLHS();
            const auto& rhs = getRHS();
            const auto numThread = std::min((unsigned int)rhs.getSize(), CUDADevAttr::DefaultThreadsPerBlock);
            const unsigned int numBlockX = 1;
            const unsigned int numBlockY = lhs.getRow();
            const unsigned int numBlockZ = lhs.getCol();
            CUDAExecutor::launch(assign_add_general, KernelConfig({numBlockX, numBlockY, numBlockZ}, numThread));
        }
    }

    template<Matrix M1, Matrix M2>
    __device__ auto device_obj<Kronecker<M1, M2>>::calc(size_t row, size_t col) const -> T {
        size_t row1 = row / getRHS().getRow();
        size_t row2 = row % getRHS().getRow();
        size_t col1 = col / getRHS().getCol();
        size_t col2 = col % getRHS().getCol();
        return getLHS().calc(row1, col1) * getRHS().calc(row2, col2);
    }

    template<Matrix M1, Matrix M2>
    __host__ __device__ size_t device_obj<Kronecker<M1, M2>>::getRow() const noexcept {
        return getLHS().getRow() * getRHS().getRow();
    }

    template<Matrix M1, Matrix M2>
    __host__ __device__ size_t device_obj<Kronecker<M1, M2>>::getCol() const noexcept {
        return getLHS().getCol() * getRHS().getCol();
    }

    template<Matrix M1, Matrix M2>
    void device_obj<Kronecker<M1, M2>>::assign_identity(Matrix auto&& target) const {
        target.zeros();
        const auto& lhs = getLHS();
        const auto& rhs = getRHS();
        for (size_t r = 0; r < lhs.getRow(); ++r) {
            size_t offsetR = r * rhs.getRow();
            rhs.assign(target.block(offsetR, rhs.getRow(), offsetR, rhs.getCol()));
        }
    }

    template<Matrix M1, Matrix M2>
    void device_obj<Kronecker<M1, M2>>::assign_add_identity(Matrix auto&& target) const {
        const auto& lhs = getLHS();
        const auto& rhs = getRHS();
        for (size_t r = 0; r < lhs.getRow(); ++r) {
            size_t offsetR = r * rhs.getRow();
            rhs.assign_add(target.block(offsetR, rhs.getRow(), offsetR, rhs.getCol()));
        }
    }

    template<Matrix M1, Matrix M2>
    void device_obj<Kronecker<M1, M2>>::assign_diagonal(Matrix auto&& target) const {
        target.zeros();
        const auto& lhs = getLHS();
        const auto& rhs = getRHS();
        for (size_t r = 0; r < lhs.getRow(); ++r) {
            size_t offsetR = r * rhs.getRow();
            (rhs * lhs.calc(r, r)).assign(target.block(offsetR, rhs.getRow(), offsetR, rhs.getCol()));
        }
    }

    template<Matrix M1, Matrix M2>
    void device_obj<Kronecker<M1, M2>>::assign_add_diagonal(Matrix auto&& target) const {
        const auto& lhs = getLHS();
        const auto& rhs = getRHS();
        for (size_t r = 0; r < lhs.getRow(); ++r) {
            size_t offsetR = r * rhs.getRow();
            (rhs * lhs.calc(r, r)).assign_add(target.block(offsetR, rhs.getRow(), offsetR, rhs.getCol()));
        }
    }

    template<Matrix M1, Matrix M2>
    __host__ __device__ auto&& device_obj<Kronecker<M1, M2>>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M1>(self.m1.getDerived());
    }

    template<Matrix M1, Matrix M2>
    __host__ __device__ auto&& device_obj<Kronecker<M1, M2>>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M2>(self.m2.getDerived());
    }

    template<Matrix M1, Matrix M2>
    [[nodiscard]] __host__ __device__ auto kronecker(M1&& m1, M2&& m2) noexcept requires(DeviceObj<M1> && DeviceObj<M2>) {
        using RetTy = device_obj<Kronecker<remove_device_obj_t<M1&&>, remove_device_obj_t<M2&&>>>;
        return RetTy(std::forward<M1>(m1), std::forward<M2>(m2));
    }
}

namespace Physica {
    template<Matrix M1, Matrix M2>
    class Traits<device_obj<Kronecker<M1, M2>>> : public Traits<Kronecker<M1, M2>> {};
}
