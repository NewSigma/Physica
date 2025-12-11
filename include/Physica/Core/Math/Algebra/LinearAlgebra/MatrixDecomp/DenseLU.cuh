/*
 * Copyright 2025 Weibo He.
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

#include "DenseLU.h"
#include "Physica/Core/Exception/CUDA/cuSolver.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"

namespace Physica {
    template<Scalar T, bool Pivot>
    class device_obj<DenseLU<T, Pivot>> {
        using host_obj = DenseLU<T, Pivot>;
        using This = device_obj<host_obj>;
        using PermType = std::conditional<Pivot, device_obj<Array<int64_t>>, PlainStruct<void>>::type;

        using Tr = T::RealType;
        constexpr static cudaDataType DataType = CUDAContext::getDataType<T>();
    private:
        device_obj<MatrixND<T>> working;
        [[no_unique_address]] PermType perm;
        device_obj<Array<std::byte>> deviceBuffer;
        Array<std::byte> hostBuffer;
        device_obj<Array<int>> err = device_obj<Array<int>>(1);
    public:
        device_obj() = default;
        device_obj(size_t size);
        device_obj(const Matrix auto& source);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void compute(const Matrix auto& source);
        void solve(device_obj<MatrixND<T>>& rhs);

        [[nodiscard]] Tr lnAbsDet() const noexcept;
        [[nodiscard]] T sgndet() const noexcept;

        void resize(size_t size);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const auto& getMatrixLU() const noexcept { return working; }
        [[nodiscard]] size_t getOrder() const noexcept { return working.getRow(); }
        [[nodiscard]] size_t getRow() const noexcept { return getOrder(); }
        [[nodiscard]] size_t getCol() const noexcept { return getOrder(); }
    private:
        [[nodiscard]] constexpr static cudaDataType getDataType() noexcept;
    };

    template<Scalar T, bool Pivot>
    device_obj<DenseLU<T, Pivot>>::device_obj(size_t size) {
        resize(size);
    }

    template<Scalar T, bool Pivot>
    device_obj<DenseLU<T, Pivot>>::device_obj(const Matrix auto& source) : device_obj(source.getRow()) {
        compute(source);
    }

    template<Scalar T, bool Pivot>
    void device_obj<DenseLU<T, Pivot>>::compute(const Matrix auto& source) {
        assert(source.isSquare());
        source.assign(working);

        auto& ctx = CUDAContext::getInstance();
        int64_t m = getOrder();
        int64_t n = m;
        constexpr cudaDataType dataTypeA = DataType;
        void* A = working.data();
        int64_t lda = m;
        int64_t* ipiv = nullptr;
        if constexpr (Pivot)
            ipiv = perm.data();
        constexpr cudaDataType computeType = DataType;
        void* bufferOnDevice = deviceBuffer.data();
        size_t workspaceInBytesOnDevice = deviceBuffer.getLength();
        void* bufferOnHost = hostBuffer.getLength() == 0 ? nullptr : hostBuffer.data();
        size_t workspaceInBytesOnHost = hostBuffer.getLength();
        int* info = err.data();
        check(cusolverDnXgetrf(ctx, ctx, m, n, dataTypeA, A, lda, ipiv, computeType, bufferOnDevice, workspaceInBytesOnDevice, bufferOnHost, workspaceInBytesOnHost, info));
    }

    template<Scalar T, bool Pivot>
    void device_obj<DenseLU<T, Pivot>>::solve(device_obj<MatrixND<T>>& rhs) {
        auto& ctx = CUDAContext::getInstance();
        constexpr auto trans = CUBLAS_OP_N;
        int64_t n = getOrder();
        int64_t nrhs = rhs.getCol();
        constexpr auto dataTypeA = DataType;
        void* A = working.data();
        int64_t lda = n;
        constexpr auto dataTypeB = DataType;
        int64_t* ipiv = nullptr;
        if constexpr (Pivot)
            ipiv = perm.data();
        void* B = rhs.data();
        int64_t ldb = n;
        int* info = err.data();
        check(cusolverDnXgetrs(ctx, ctx, trans, n, nrhs, dataTypeA, A, lda, ipiv, dataTypeB, B, ldb, info));
    }

    template<Scalar T, bool Pivot>
    auto device_obj<DenseLU<T, Pivot>>::lnAbsDet() const noexcept -> Tr {
        return ln(abs(working.diag())).sum();
    }

    template<Scalar T, bool Pivot>
    T device_obj<DenseLU<T, Pivot>>::sgndet() const noexcept {
        return unit(working.diag()).prod();
    }

    template<Scalar T, bool Pivot>
    void device_obj<DenseLU<T, Pivot>>::resize(size_t size) {
        working.resize(size, size);
        if constexpr (Pivot)
            perm.resize(size);
        size_t deviceBufferSize = 0, hostBufferSize = 0;
        {
            auto& ctx = CUDAContext::getInstance();
            int64_t m = getOrder();
            int64_t n = m;
            constexpr cudaDataType dataTypeA = DataType;
            const void* A = working.data();
            int64_t lda = m;
            constexpr cudaDataType computeType = DataType;
            check(cusolverDnXgetrf_bufferSize(ctx, ctx, m, n, dataTypeA, A, lda, computeType, &deviceBufferSize, &hostBufferSize));
        }
        deviceBuffer.resize(deviceBufferSize);
        hostBuffer.resize(hostBufferSize);
    }

    template<Scalar T, bool Pivot>
    void device_obj<DenseLU<T, Pivot>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        working.swap(obj.working);
        perm.swap(obj.perm);
        deviceBuffer.swap(obj.deviceBuffer);
        hostBuffer.swap(obj.hostBuffer);
        err.swap(obj.err);
    }
}
