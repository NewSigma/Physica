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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/PermMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/SparseMatrix.h"

namespace Physica {
    template<Scalar T>
    class SparseLU {
        using This = SparseLU<T>;
        using Tr = T::RealType;
        enum Phase : int64_t {
            Analyse = 11,
            Factorize = 22,
            AnalyseFactorize = 12,
            AnalyseFactorizeSolve = 13,
            Solve = 33,
            Destroy = -1,
        };

        SparseMatrix<T> spmat;

        VectorND<T> diags;
        VectorND<T> buffer;

        Array<void*, 64> handles;
        Array<int64_t, 64> params;
    public:
        SparseLU();
        SparseLU(size_t size, bool needDiag);
        SparseLU(SparseMatrix<T> source, bool needDiag);
        SparseLU(const This&) = delete;
        SparseLU(This&& obj) noexcept;
        ~SparseLU();
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void compute(const Matrix auto& source);
        void analyse(const Matrix auto& source);
        void factorize(const Matrix auto& source);
        [[nodiscard]] MatrixND<T> solve(const MatrixND<T>& rhs);
        [[nodiscard]] MatrixND<T> solve_mkl(const MatrixND<T>& rhs);

        [[nodiscard]] Tr lnAbsDet() const noexcept;
        [[nodiscard]] T sgndet() const noexcept;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getOrder() const noexcept { return spmat.getRow(); }
        [[nodiscard]] bool getNeedDiag() const noexcept { return params[55]; }
    private:
        void initialize(bool needDiag) noexcept;
        void runPardiso(Phase p, int64_t nrhs = 0, T* b = nullptr, T* x = nullptr);
        /* Static members */
        constexpr static bool isFactorizePhase(Phase p) noexcept;
        constexpr static bool isSolvePhase(Phase p) noexcept;
        constexpr static int getMatrixType() noexcept;
    };

    template<Scalar T>
    SparseLU<T>::SparseLU() {
        handles.zeros();
        params.zeros();
    }

    template<Scalar T>
    SparseLU<T>::SparseLU(This&& obj) noexcept
            : spmat(std::move(obj.spmat))
            , diags(std::move(obj.diags))
            , buffer(std::move(obj.buffer))
            , handles(std::move(obj.handles))
            , params(std::move(obj.params)) {
        handles.zeros();
        params.zeros();
    }

    template<Scalar T>
    SparseLU<T>::SparseLU(size_t size, bool needDiag) : spmat(size, size) {
        if (needDiag) {
            diags.resize(getOrder());
            buffer.resize(getOrder());
        }
        initialize(needDiag);
    }

    template<Scalar T>
    SparseLU<T>::SparseLU(SparseMatrix<T> source, bool needDiag) : SparseLU(source.getRow(), needDiag) {
        assert(source.isSquare());
        spmat = std::move(source);
        runPardiso(Phase::AnalyseFactorize);
    }

    template<Scalar T>
    SparseLU<T>::~SparseLU() {
        runPardiso(Phase::Destroy);
    }

    template<Scalar T>
    MatrixND<T> SparseLU<T>::solve(const MatrixND<T>& rhs) {
        if constexpr (HasMKL())
            return solve_mkl(rhs);
        else
            noImpl(__func__);
    }

    template<Scalar T>
    auto SparseLU<T>::lnAbsDet() const noexcept -> Tr {
        assert(getNeedDiag());
        return ln(abs(diags)).sum();
    }

    template<Scalar T>
    T SparseLU<T>::sgndet() const noexcept {
        assert(getNeedDiag());
        noImpl("[Error]: MKL does not provide API for sgndet implementation");
    }

    template<Scalar T>
    void SparseLU<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        spmat.swap(obj.spmat);
        diags.swap(obj.diags);
        buffer.swap(obj.buffer);
        handles.swap(obj.handles);
        params.swap(obj.params);
    }
}

#ifdef PHYSICA_MKL
    #include "SparseLUImpl/SparseLU_MKL.h"
#endif
