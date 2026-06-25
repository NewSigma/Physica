/*
 * Copyright 2026 Weibo He.
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

#include "MatrixTrig.cuh"
#include "Physica/Core/Utils/Container/Array.cuh"

namespace Physica {
    template<Matrix M> requires(instanceof_tx<M, MatrixTrig>)
    class device_obj<Inverse<M>> : public device_obj<RValueMatrix<Inverse<M>>> {
        using host_obj = Inverse<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using Ref = add_device_obj<M>::type;
    protected:
        using typename Base::T;
    private:
        constexpr static cudaDataType DataType = CUDAContext::getDataType<T>();

        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> trig;

        device_obj<Array<std::byte>> deviceBuffer;
        Array<std::byte> hostBuffer;
        device_obj<Array<int>> err = device_obj<Array<int>>(1);
    public:
        explicit device_obj(Ref trig);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Operations */
        void assign(Matrix auto& target) const;
        void assign_cusolver(Matrix auto& target) const;
        /* Getters */
        [[nodiscard]] auto&& getExpr(this auto&&) noexcept;
        [[nodiscard]] size_t getOrder() const noexcept { return getExpr().getOrder(); }
    private:
        void allocate();
    };

    template<Matrix M> requires(instanceof_tx<M, MatrixTrig>)
    device_obj<Inverse<M>>::device_obj(Ref trig) : trig(asStruct(trig)) {
        allocate();
    }

    template<Matrix M> requires(instanceof_tx<M, MatrixTrig>)
    void device_obj<Inverse<M>>::assign(Matrix auto& target) const {
        using Expr = Traits<M>::ExprType;
        if constexpr (Internal::EnableLAPACK<Expr, decltype(target)>::value)
            assign_cusolver(target);
        else
            noImpl();
    }

    template<Matrix M> requires(instanceof_tx<M, MatrixTrig>)
    void device_obj<Inverse<M>>::assign_cusolver(Matrix auto& target) const {
        using Tm = decltype(std::declval<T>().toCUDA());
        getExpr().assign(target);

        auto& ctx = CUDAContext::getInstance();
        constexpr auto Uplo = Traits<M>::Upper ? cublasFillMode_t::CUBLAS_FILL_MODE_UPPER : cublasFillMode_t::CUBLAS_FILL_MODE_LOWER;
        constexpr auto Diag = Traits<M>::Unit ? cublasDiagType_t::CUBLAS_DIAG_UNIT : cublasDiagType_t::CUBLAS_DIAG_NON_UNIT;
        size_t n = getOrder();
        auto* A = reinterpret_cast<Tm*>(target.data());
        void* bufferOnDevice = deviceBuffer.data();
        size_t workspaceInBytesOnDevice = deviceBuffer.getLength();
        void* bufferOnHost = hostBuffer.getLength() == 0 ? nullptr : hostBuffer.data();
        size_t workspaceInBytesOnHost = hostBuffer.getLength();
        int* info = err.data();
        check(cusolverDnXtrtri(ctx, Uplo, Diag, n, DataType, A, n, bufferOnDevice, workspaceInBytesOnDevice, bufferOnHost, workspaceInBytesOnHost, info));
    }

    template<Matrix M> requires(instanceof_tx<M, MatrixTrig>)
    auto&& device_obj<Inverse<M>>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref>(self.trig.getDerived());
    }

    template<Matrix M> requires(instanceof_tx<M, MatrixTrig>)
    void device_obj<Inverse<M>>::allocate() {
        size_t workspaceInBytesOnDevice{}, workspaceInBytesOnHost{};
        {
            auto& ctx = CUDAContext::getInstance();
            constexpr auto Uplo = Traits<M>::Upper ? cublasFillMode_t::CUBLAS_FILL_MODE_UPPER : cublasFillMode_t::CUBLAS_FILL_MODE_LOWER;
            constexpr auto Diag = Traits<M>::Unit ? cublasDiagType_t::CUBLAS_DIAG_UNIT : cublasDiagType_t::CUBLAS_DIAG_NON_UNIT;
            size_t n = getOrder();
            check(cusolverDnXtrtri_bufferSize(ctx, Uplo, Diag, n, DataType, nullptr, n, &workspaceInBytesOnDevice, &workspaceInBytesOnHost));
        }
        deviceBuffer.resize(workspaceInBytesOnDevice);
        hostBuffer.resize(workspaceInBytesOnHost);
    }
}
