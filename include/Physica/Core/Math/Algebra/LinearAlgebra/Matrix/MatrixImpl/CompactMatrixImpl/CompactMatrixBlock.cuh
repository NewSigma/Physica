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

#include "CompactMatrixBlock.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/CompactVector.cuh"
#include "Physica/Core/Utils/Empty.h"

namespace Physica {
    template<Matrix M, size_t Row, size_t Col>
    class device_obj<CompactMatrixBlock<M, Row, Col>> : public device_obj<CompactVector<CompactMatrixBlock<M, Row, Col>>> {
        using host_obj = CompactMatrixBlock<M, Row, Col>;
        using This = device_obj<host_obj>;
        using Base = device_obj<CompactVector<host_obj>>;
        using Ref = add_device_obj<M>::type;
        using MaybeRowCount = std::conditional_t<Row == Dynamic, size_t, Empty>;
        using MaybeColCount = std::conditional_t<Col == Dynamic, size_t, Empty>;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
        size_t fromRow;
        [[no_unique_address]] MaybeRowCount rowCount;
        size_t fromCol;
        [[no_unique_address]] MaybeColCount colCount;
    public:
        __host__ __device__ device_obj(Ref mat, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount);
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This& m);
        This& operator=(This&& m) noexcept;
        using Base::operator=;
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize(size_t length);
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr size_t getLength() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle(this auto&& self) noexcept;
    };

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ device_obj<CompactMatrixBlock<M, Row, Col>>::device_obj(Ref mat, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount)
            : mat(asStruct(mat))
            , fromRow(fromRow)
            , rowCount(rowCount)
            , fromCol(fromCol)
            , colCount(colCount) {
        assert(fromRow < mat.getRow());
        assert(fromCol < mat.getCol());
        assert(fromRow + rowCount <= mat.getRow());
        assert(fromCol + colCount <= mat.getCol());
        if constexpr (Row != Dynamic)
            assert(Row == rowCount);
        if constexpr (Col != Dynamic)
            assert(Col == colCount);
    }

    template<Matrix M, size_t Row, size_t Col>
    auto device_obj<CompactMatrixBlock<M, Row, Col>>::operator=(const This& m) -> This& {
        Base::operator=(m);
        return *this;
    }

    template<Matrix M, size_t Row, size_t Col>
    auto device_obj<CompactMatrixBlock<M, Row, Col>>::operator=(This&& m) noexcept -> This& {
        return *this = m;
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ void device_obj<CompactMatrixBlock<M, Row, Col>>::resize([[maybe_unused]] size_t length) {
        assert(length == getLength());
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ constexpr size_t device_obj<CompactMatrixBlock<M, Row, Col>>::getLength() const noexcept {
        if constexpr (Row == 1) {
            if constexpr (Col == Dynamic)
                return colCount;
            else
                return Col;
        }
        else {
            if constexpr (Row == Dynamic)
                return rowCount;
            else
                return Row;
        }
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrixBlock<M, Row, Col>>::data_handle(this auto&& self) noexcept {
        return self.mat.getDerived().data_ptr(self.fromRow, self.fromCol);
    }
}

namespace Physica {
    template<Matrix M, size_t Row, size_t Col>
    class Traits<device_obj<CompactMatrixBlock<M, Row, Col>>> : public Traits<CompactMatrixBlock<M, Row, Col>> {
        static_assert(!is_device_obj<M>::value);
    };
}
