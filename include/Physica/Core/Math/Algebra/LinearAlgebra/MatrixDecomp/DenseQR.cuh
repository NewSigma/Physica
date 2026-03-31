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

#include "DenseQR.h"
#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"
#include "Physica/Core/Utils/Container/Array.cuh"

namespace Physica {
    template<Scalar T>
    class device_obj<DenseQR<T>> {
        static_assert(T::Prec == Float32 || T::Prec == Float64);
        using host_obj = DenseQR<T>;
        using This = device_obj<host_obj>;
        using DeviceVector = device_obj<VectorND<T>>;
        using Tc = T::ComplexType;
        constexpr static bool isComplex = T::isComplex();
    private:
        device_obj<MatrixND<T>> working;
        DeviceVector taus;
        DeviceVector vecD;
        device_obj<Array<std::byte>> deviceBuffer;
        Array<std::byte> hostBuffer;
    public:
        device_obj() = default;
        device_obj(size_t row, size_t col);
        device_obj(const Matrix auto& source);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void compute(const Matrix auto& source);

        [[nodiscard]] __device__ T calcDetQ() const;
        void toQDT();

        void resize(size_t row, size_t col);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ auto&& getWorking(this auto&& self) noexcept;
        [[nodiscard]] __host__ __device__ const auto& getTaus() const noexcept { return taus; }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return working.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return working.getCol(); }
        [[nodiscard]] device_obj<MatrixND<T>> getMatrixQ();
        [[nodiscard]] __host__ __device__ auto getMatrixR() const noexcept;
    };

    template<Scalar T>
    device_obj<DenseQR<T>>::device_obj(size_t row, size_t col) {
        resize(row, col);
    }

    template<Scalar T>
    device_obj<DenseQR<T>>::device_obj(const Matrix auto& source) : device_obj(source.getRow(), source.getCol()) {
        compute(source);
    }

    template<Scalar T>
    void device_obj<DenseQR<T>>::compute(const Matrix auto& source) {
        source.assign(working);

        auto& ctx = CUDAContext::getInstance();
        const size_t m = getRow();
        const size_t n = getCol();
        constexpr cudaDataType dataType = CUDAContext::getDataType<T>();
        void* A = working.data();
        const size_t lda = m;
        void* tau = taus.data();
        int info;
        void* bufferOnHost = hostBuffer.getLength() == 0 ? nullptr : hostBuffer.data();
        check(cusolverDnXgeqrf(ctx, ctx, m, n, dataType, A, lda, dataType, tau, dataType, deviceBuffer.data(), deviceBuffer.getLength(), bufferOnHost, hostBuffer.getLength(), &info));
    }

    template<Scalar T>
    __device__ T device_obj<DenseQR<T>>::calcDetQ() const {
        int sign = 0;
        for (auto tau : taus)
            sign += tau.isPositive();
        return T(sign % 2 == 0 ? 1.0 : -1.0);
    }

    template<Scalar T>
    void device_obj<DenseQR<T>>::toQDT() {
        const size_t length = taus.getLength();
        vecD.resize(length);

        const int numThread = std::min<size_t>(length, CUDADevAttr::DefaultThreadsPerBlock);
        CUDAExecutor::launch([working_ = asStruct(working), vecD_ = asStruct(vecD)] __device__() mutable {
            auto& working = working_.getDerived();
            auto& vecD = vecD_.getDerived();
            const int row = blockIdx.x;
            if (working[row, row].isZero()) {
                if (isZeroThread())
                    vecD[row] = 1;
                return;
            }

            __shared__ T factor;
            if (isZeroThread()) {
                vecD[row] = working[row, row];
                factor = reciprocal(vecD[row]);
            }
            __syncthreads();

            for (int i = threadIdx.x + row; i < vecD.getLength(); i += blockDim.x)
                working[row, i] *= factor;
        }, KernelConfig(length, numThread));
    }

    template<Scalar T>
    void device_obj<DenseQR<T>>::resize(size_t row, size_t col) {
        working.resize(row, col);
        taus.resize(std::min(row, col));
        taus.zeros(); // For historic reason, BLAS-like interface will allocate a unused element

        auto& ctx = CUDAContext::getInstance();
        size_t deviceBufferSize = 0, hostBufferSize = 0;
        /* Buffer size check 1 */ {
            const size_t m = getRow();
            const size_t n = getCol();
            constexpr cudaDataType dataType = CUDAContext::getDataType<T>();
            const size_t lda = m;
            check(cusolverDnXgeqrf_bufferSize(ctx, ctx, m, n, dataType, nullptr, lda, dataType, nullptr, dataType, &deviceBufferSize, &hostBufferSize));
        }
        /* Buffer size check 2 */ {
            auto& ctx = CUDAContext::getInstance();
            const size_t m = getRow();
            const size_t n = m;
            const size_t k = getCol();
            const size_t lda = m;
            int lwork;
            if constexpr (T::Prec == Float32)
                check(cusolverDnSorgqr_bufferSize(ctx, m, n, k, nullptr, lda, nullptr, &lwork));
            else
                check(cusolverDnDorgqr_bufferSize(ctx, m, n, k, nullptr, lda, nullptr, &lwork));
            deviceBufferSize = std::max(deviceBufferSize, (size_t)lwork);
        }
        deviceBuffer.resize(deviceBufferSize);
        hostBuffer.resize(hostBufferSize);
    }

    template<Scalar T>
    void device_obj<DenseQR<T>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        working.swap(obj.working);
        taus.swap(obj.taus);
        hostBuffer.swap(obj.hostBuffer);
        deviceBuffer.swap(obj.deviceBuffer);
    }

    template<Scalar T>
    __host__ __device__ auto&& device_obj<DenseQR<T>>::getWorking(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).working;
    }

    template<Scalar T>
    auto device_obj<DenseQR<T>>::getMatrixQ() -> device_obj<MatrixND<T>> {
        static_assert(!isComplex, "[Error]: Not implemented");
        using Tm = decltype(std::declval<T>().toCUDA());

        auto& ctx = CUDAContext::getInstance();
        const size_t m = getRow();
        const size_t n = m;
        const size_t k = getCol();
        const size_t lda = m;
        device_obj<MatrixND<T>> result(m, lda);
        result.leftCols(k) = working;
        auto* A = reinterpret_cast<Tm*>(result.data());
        auto* tau = reinterpret_cast<const Tm*>(taus.data());
        auto* work = reinterpret_cast<Tm*>(deviceBuffer.data());
        size_t lwork = deviceBuffer.getLength();

        int info;
        if constexpr (T::Prec == Float32)
            check(cusolverDnSorgqr(ctx, m, n, k, A, lda, tau, work, lwork, &info));
        else
            check(cusolverDnDorgqr(ctx, m, n, k, A, lda, tau, work, lwork, &info));
        return result;
    }

    template<Scalar T>
    __host__ __device__ auto device_obj<DenseQR<T>>::getMatrixR() const noexcept {
        return working.triu();
    }
}
