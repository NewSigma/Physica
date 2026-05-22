/*
 * Copyright 2023-2026 Weibo He.
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

#include "LValueMatrix.h"
#include "CompactMatrixImpl/CompactMatrixBlock.h"

namespace Physica {
    template<class> class FlattenC;
    /**
     * \class CompactMatrix has its elements on major direction distributed Compactly.
     */
    template<class Derived>
    class CompactMatrix : public LValueMatrix<Derived> {
        using Base = LValueMatrix<Derived>;
        using This = CompactMatrix<Derived>;
    public:
        using typename Base::ScalarType;
        using Base::isComplex;
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
    public:
        ~CompactMatrix() = default;
        /* Operators */
        This& operator=(const This& obj) = delete;
        This& operator=(This&& obj) noexcept = delete;
        using Base::operator=;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Matrix auto&& m) const noexcept;
        void assign_mkl(Matrix auto&& target) const noexcept;
        template<ExecutePolicy P = Sequential>
        void assign_base(Matrix auto&& __restrict target) const __restrict noexcept;

        template<Matrix M> void toDevice(device_obj<CompactMatrix<M>>& obj) const;
        template<Matrix M> void toDeviceAsync(device_obj<CompactMatrix<M>>& obj) const;
        [[nodiscard]] auto toNumpy() const;

        [[nodiscard]] auto row(this auto&&, size_t r) noexcept;
        [[nodiscard]] auto col(this auto&&, size_t c) noexcept;
        template<size_t Row = Dynamic>
        [[nodiscard]] auto rows(this auto&&, size_t fromRow, size_t rowCount) noexcept;
        template<size_t Row = Dynamic>
        [[nodiscard]] auto topRows(this auto&&, size_t to) noexcept;
        template<size_t Row = Dynamic>
        [[nodiscard]] auto bottomRows(this auto&&, size_t from) noexcept;
        template<size_t Col = Dynamic>
        [[nodiscard]] auto cols(this auto&&, size_t fromCol, size_t colCount) noexcept;
        template<size_t Col = Dynamic>
        [[nodiscard]] auto leftCols(this auto&&, size_t to) noexcept;
        template<size_t Col = Dynamic>
        [[nodiscard]] auto rightCols(this auto&&, size_t from) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto topLeftCorner(this auto&&, size_t toRow, size_t toCol) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto topLeftCorner(this auto&&, size_t to) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto topRightCorner(this auto&&, size_t toRow, size_t fromCol) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto bottomLeftCorner(this auto&&, size_t fromRow, size_t toCol) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto bottomRightCorner(this auto&&, size_t fromRow, size_t fromCol) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto bottomRightCorner(this auto&&, size_t from) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto block(this auto&&, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept;

        [[nodiscard]] auto flatten(this auto&&) noexcept;

        void zeros() noexcept;
        void read(const auto& obj) noexcept;
        template<RNG R>
        void random_uniform();
        template<RNG R>
        void random_normal();
        [[nodiscard]] VectorND<T> balance_mkl();
        const H5DataSet<2> read(const H5Loc& loc, const char* name);
        H5DataSet<2> write(H5Loc& loc, const char* name) const;
        /* Getters */
        [[nodiscard]] auto data() noexcept;
        [[nodiscard]] auto data() const noexcept;
        [[nodiscard]] auto data_ptr(this auto&&, size_t r, size_t c) noexcept;
    protected:
        CompactMatrix() = default;
        CompactMatrix(const This&) = default;
        CompactMatrix(This&&) noexcept = default;
    };
}

namespace Physica {
    template<class Derived>
    class Traits<CompactMatrix<Derived>> : public Traits<Derived> {};
}

#include "CompactMatrixImpl/CompactMatrixImpl.h"
#ifdef PHYSICA_MKL
    #include "CompactMatrixImpl/CompactMatrix_MKL.h"
#endif
#include "CompactMatrixImpl/Flatten.h"
#include "CompactMatrixImpl/Transpose.h"
#include "CompactMatrixImpl/ReshapedVector.h"
