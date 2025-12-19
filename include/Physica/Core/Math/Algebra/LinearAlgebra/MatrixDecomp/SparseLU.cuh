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

#include "SparseLU.h"
#include "Physica/Core/Parallel/CUDAContext.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/SparseMatrix.cuh"

namespace Physica {
    template<Scalar T>
    class device_obj<SparseLU<T>> {
        using host_obj = SparseLU<T>;
        using This = device_obj<host_obj>;
        struct ConfigDeleter {
            void operator()(cudssConfig_t ptr) const { cudssConfigDestroy(ptr); } 
        };

        struct DataDeleter {
            void operator()(cudssData_t ptr) const { cudssDataDestroy(CUDAContext::getInstance(), ptr); } 
        };

        struct MatrixDeleter {
            void operator()(cudssMatrix_t ptr) const { cudssMatrixDestroy(ptr); }
        };

        using MatrixCUDSS = std::unique_ptr<cudssMatrix, MatrixDeleter>;
        using Tr = T::RealType;
        using Trv = Tr::ValueType;
    private:
        device_obj<SparseMatrix<T>> spmat;
        device_obj<VectorND<T>> diag;
        device_obj<Array<size_t>> permRow;
        device_obj<Array<size_t>> permCol;

        std::unique_ptr<cudssConfig, ConfigDeleter> config = nullptr;
        std::unique_ptr<cudssData, DataDeleter> data = nullptr;
        MatrixCUDSS sparseA = nullptr;
    public:
        device_obj() = default;
        device_obj(size_t size);
        device_obj(const Matrix auto& source);
        device_obj(const This&) = delete;
        device_obj(This&& obj) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void compute(const Matrix auto& source);
        void refactorize(const SparseMatrix<T>& source);
        void solve(device_obj<MatrixND<T>>& rhs);

        [[nodiscard]] Tr lnAbsDet();
        [[nodiscard]] T sgndet();

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getOrder() const noexcept { return spmat.getRow(); }
    private:
        void initialize();
        [[nodiscard]] MatrixCUDSS makeDenseCUDSS(const device_obj<MatrixND<T>>& m);
        [[nodiscard]] MatrixCUDSS makeSparseCUDSS();
        void compute_cudss(cudssPhase_t phase, cudssMatrix_t sol = nullptr, cudssMatrix_t rhs = nullptr);
        void collectDiag();
    };

    template<Scalar T>
    device_obj<SparseLU<T>>::device_obj(size_t size)
            : spmat(size, size)
            , diag(size)
            , permRow(size)
            , permCol(size) {
        initialize();
    }

    template<Scalar T>
    device_obj<SparseLU<T>>::device_obj(const Matrix auto& source) : device_obj(source.getRow()) {
        compute(source);
    }

    template<Scalar T>
    void device_obj<SparseLU<T>>::compute(const Matrix auto& source) {
        assert(source.isSquare());
        using M = std::remove_cvref<decltype(source)>::type;
        if constexpr (std::same_as<M, SparseMatrix<T>>)
            source.toDeviceAsync(spmat);
        else if constexpr (CUDA<M>)
            spmat = source;
        else
            spmat = SparseMatrix<T>(source);
        sparseA = makeSparseCUDSS();
        compute_cudss(cudssPhase_t(CUDSS_PHASE_ANALYSIS | CUDSS_PHASE_FACTORIZATION));
        collectDiag();
    }

    template<Scalar T>
    void device_obj<SparseLU<T>>::refactorize(const SparseMatrix<T>& source) {
        assert(getOrder() == source.getRow() && "[Error]: Matrix size do not match");
        assert(source.isSquare());
        source.getElements().toDeviceAsync(spmat.getElements());
        compute_cudss(CUDSS_PHASE_REFACTORIZATION);
        collectDiag();
    }

    template<Scalar T>
    void device_obj<SparseLU<T>>::solve(device_obj<MatrixND<T>>& rhs) {
        auto rhsMatrix = makeDenseCUDSS(rhs);
        compute_cudss(CUDSS_PHASE_SOLVE, rhsMatrix.get(), rhsMatrix.get());
    }

    template<Scalar T>
    auto device_obj<SparseLU<T>>::lnAbsDet() -> Tr {
        return ln(abs(diag)).sum();
    }

    template<Scalar T>
    T device_obj<SparseLU<T>>::sgndet() {
        Trv permDet = PermMatrix<Trv>(permRow.toHost()).det() * PermMatrix<Trv>(permCol.toHost()).det();
        return unit(diag).prod() * permDet;
    }

    template<Scalar T>
    void device_obj<SparseLU<T>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        spmat.swap(obj.spmat);
        diag.swap(obj.diag);
        permRow.swap(obj.permRow);
        permCol.swap(obj.permCol);
        config.swap(obj.config);
        data.swap(obj.data);
        sparseA.swap(obj.sparseA);
    }

    template<Scalar T>
    void device_obj<SparseLU<T>>::initialize() {
        {
            cudssConfig_t pConfig = nullptr;
            check(cudssConfigCreate(&pConfig));
            config.reset(pConfig);
        }
        {
            auto& ctx = CUDAContext::getInstance();
            cudssData_t pData = nullptr;
            check(cudssDataCreate(ctx, &pData));
            data.reset(pData);
        }
        cudssAlgType_t value = CUDSS_ALG_1;
        check(cudssConfigSet(config.get(), CUDSS_CONFIG_REORDERING_ALG, &value, sizeof(value)));
    }

    template<Scalar T>
    auto device_obj<SparseLU<T>>::makeDenseCUDSS(const device_obj<MatrixND<T>>& m) -> MatrixCUDSS {
        cudssMatrix_t handle = nullptr;
        int64_t nrows = m.getRow();
        int64_t ncols = m.getCol();
        int64_t ld = m.getMaxMinor();
        auto* values = m.data();
        constexpr auto valueType = CUDAContext::getDataType<T>();
        constexpr auto layout = cudssLayout_t::CUDSS_LAYOUT_COL_MAJOR;
        check(cudssMatrixCreateDn(&handle, nrows, ncols, ld, (void*)values, valueType, layout));
        return MatrixCUDSS(handle);
    }

    template<Scalar T>
    auto device_obj<SparseLU<T>>::makeSparseCUDSS() -> MatrixCUDSS {
        cudssMatrix_t handle = nullptr;
        int64_t n = getOrder();
        int64_t nnz = spmat.getNumNonZero();
        auto* rowStart = spmat.getMajorStarts().data();
        auto* colIndices = spmat.getMinorIndexes().data();
        auto* values = spmat.getElements().data();
        constexpr auto indexType = cudaDataType::CUDA_R_64I;
        constexpr auto valueType = CUDAContext::getDataType<T>();
        constexpr auto mtype = cudssMatrixType_t::CUDSS_MTYPE_GENERAL;
        constexpr auto mview = cudssMatrixViewType_t::CUDSS_MVIEW_FULL;
        constexpr auto indexBase = cudssIndexBase_t::CUDSS_BASE_ZERO;
        check(cudssMatrixCreateCsr(&handle, n, n, nnz, (void*)rowStart, nullptr, (void*)colIndices, (void*)values, indexType, valueType, mtype, mview, indexBase));
        return MatrixCUDSS(handle);
    }

    template<Scalar T>
    void device_obj<SparseLU<T>>::compute_cudss(cudssPhase_t phase, cudssMatrix_t sol, cudssMatrix_t rhs) {
        auto& ctx = CUDAContext::getInstance();
        check(cudssExecute(ctx, phase, config.get(), data.get(), sparseA.get(), sol, rhs));
    }

    template<Scalar T>
    void device_obj<SparseLU<T>>::collectDiag() {
        auto& ctx = CUDAContext::getInstance();
        const size_t length = diag.getLength();
        {
            const size_t sizeInBytes = length * sizeof(T);
            size_t sizeWritten = 0;
            check(cudssDataGet(ctx, data.get(), cudssDataParam_t::CUDSS_DATA_DIAG, diag.data(), sizeInBytes, &sizeWritten));
            assert(sizeInBytes == sizeWritten);
        }

        const size_t sizeInBytes = length * sizeof(size_t);
        {
            size_t sizeWritten = 0;
            check(cudssDataGet(ctx, data.get(), cudssDataParam_t::CUDSS_DATA_PERM_ROW, permRow.data(), sizeInBytes, &sizeWritten));
            assert(sizeInBytes == sizeWritten);
        }
        {
            size_t sizeWritten = 0;
            check(cudssDataGet(ctx, data.get(), cudssDataParam_t::CUDSS_DATA_PERM_COL, permCol.data(), sizeInBytes, &sizeWritten));
            assert(sizeInBytes == sizeWritten);
        }
    }
}
