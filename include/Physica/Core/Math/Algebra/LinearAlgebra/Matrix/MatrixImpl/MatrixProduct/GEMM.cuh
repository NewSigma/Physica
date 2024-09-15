/*
 * Copyright 2024 Weibo He.
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

#include <Physica/PlainStruct.h>
#include <Physica/Core/Exception/NoImplException.h>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixOption.h>
#include <Physica/Core/Exception/CUDA/cuBLAS.cuh>

namespace Physica::Core {
    template<class MatrixType1, class MatrixType2>
    class device_obj<MatrixProduct<MatrixType1, MatrixType2>> : public device_obj<RValueMatrix<MatrixProduct<MatrixType1, MatrixType2>>> {
        using host_obj = MatrixProduct<MatrixType1, MatrixType2>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
    public:
        using Base::isReverseDiff;
        using typename Base::ScalarType;
        using DefaultType = typename host_obj::DefaultType::device_obj_type;
    private:
        Physica::PlainStruct<const device_obj<MatrixType1>> mat1;
        Physica::PlainStruct<const device_obj<MatrixType2>> mat2;
    public:
        __host__ __device__ device_obj(const device_obj<RValueMatrix<MatrixType1>>& mat1_, const device_obj<RValueMatrix<MatrixType2>>& mat2_);
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<class OtherDerived>
        void assignTo(device_obj<ContinuousMatrix<OtherDerived>>& target) const;
        [[nodiscard]] DefaultType compute() const { return DefaultType(*this); }
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t, size_t) const { noImpl(); }
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getRow() const { return getLHS().getRow(); }
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getColumn() const { return getRHS().getColumn(); }
        [[nodiscard]] __host__ __device__ const device_obj<MatrixType1>& getLHS() const noexcept { return mat1.getDerived(); }
        [[nodiscard]] __host__ __device__ const device_obj<MatrixType2>& getRHS() const noexcept { return mat2.getDerived(); }
    };

    template<class MatrixType1, class MatrixType2>
    __host__ __device__ device_obj<MatrixProduct<MatrixType1, MatrixType2>>::device_obj(
            const device_obj<RValueMatrix<MatrixType1>>& mat1_, const device_obj<RValueMatrix<MatrixType2>>& mat2_)
            : mat1(asStruct(mat1_.getDerived())), mat2(asStruct(mat2_.getDerived())) {
        assert(mat1_.getColumn() == mat2_.getRow());
    }

    template<class MatrixType1, class MatrixType2>
    template<class OtherDerived>
    void device_obj<MatrixProduct<MatrixType1, MatrixType2>>::assignTo(device_obj<ContinuousMatrix<OtherDerived>>& target) const {
        static_assert(MatrixOption::isColumnMatrix<MatrixType1>() && MatrixOption::isColumnMatrix<MatrixType2>(), "[Error]: cuBLAS uses column major");
        static_assert(MatrixOption::isElementMatrix<MatrixType1>() && MatrixOption::isElementMatrix<MatrixType2>(), "[Error]: cuBLAS need element storage");
        constexpr bool IsDeviceMatrix = Traits<MatrixType1>::SizeAtCompile == Dynamic && Traits<MatrixType2>::SizeAtCompile == Dynamic;
        static_assert(IsDeviceMatrix, "[Error]: Fixed matrix is on host, pass it to device before calling cuBLAS");

        using T = typename ScalarType::TrivialType;
        const auto& lhs = getLHS();
        const auto& rhs = getRHS();
        const size_t r = getRow();
        const size_t c = getColumn();
        const size_t k = lhs.getColumn();
        const ScalarType alpha = 1;
        const ScalarType beta = 0;
        auto& ctx = CUDAContext::getInstance();
        ctx.setPointerMode(false);
        if constexpr (ScalarType::Option == Float16) {
            check(cublasHgemm_64(ctx, CUBLAS_OP_N, CUBLAS_OP_N, r, c, k
                    , (T*)&alpha, (T*)lhs.data(), r, (T*)rhs.data(), k, (T*)&beta, (T*)target.data(), r));
        }
        else if constexpr (ScalarType::Option == Float32) {
            check(cublasSgemm_64(ctx, CUBLAS_OP_N, CUBLAS_OP_N, r, c, k
                    , (T*)&alpha, (T*)lhs.data(), r, (T*)rhs.data(), k, (T*)&beta, (T*)target.data(), r));
        }
        else if constexpr (ScalarType::Option == Float64) {
            check(cublasDgemm_64(ctx, CUBLAS_OP_N, CUBLAS_OP_N, r, c, k
                    , (T*)&alpha, (T*)lhs.data(), r, (T*)rhs.data(), k, (T*)&beta, (T*)target.data(), r));
        }
        else
            static_assert(ScalarType::Option == Float16, "[Error]: Unknown ScalarType");
    }

    template<class MatrixType1, class MatrixType2>
    [[nodiscard]] inline device_obj<MatrixProduct<MatrixType1, MatrixType2>>
    operator*(const device_obj<RValueMatrix<MatrixType1>>& mat1, const device_obj<RValueMatrix<MatrixType2>>& mat2) noexcept {
        return {mat1, mat2};
    }
}

namespace Physica {
    template<class MatrixType1, class MatrixType2>
    class Traits<Core::device_obj<MatrixProduct<MatrixType1, MatrixType2>>>
            : public Traits<MatrixProduct<MatrixType1, MatrixType2>> {};
}
