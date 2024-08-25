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

namespace Physica::Core {
    template<class MatrixType, class VectorType>
    class device_obj<MatrixVectorProduct<MatrixType, VectorType>>
            : public device_obj<RValueVector<MatrixVectorProduct<MatrixType, VectorType>>> {
        static_assert(MatrixType::ColumnAtCompile == VectorType::SizeAtCompile ||
                      MatrixType::ColumnAtCompile == Dynamic ||
                      VectorType::SizeAtCompile == Dynamic,
                      "[Error]: Row and column do not match in matrix-vector product");
    public:
        using host_obj = MatrixVectorProduct<MatrixType, VectorType>;
        using Base = device_obj<RValueVector<host_obj>>;
        using typename Base::ScalarType;
    private:
        using This = device_obj<host_obj>;
        Physica::PlainStruct<const device_obj<MatrixType>> mat;
        Physica::PlainStruct<const device_obj<VectorType>> vec;
    public:
        __host__ __device__ device_obj(const device_obj<RValueMatrix<MatrixType>>& mat_, const device_obj<RValueVector<VectorType>>& vec_)
                : mat(asStruct(mat_.getDerived())), vec(asStruct(vec_.getDerived())) {
            assert(mat.getDerived().getColumn() == vec.getDerived().getLength());
        }
        /* Operations */
        template<class OtherDerived>
        __device__ void assignTo(device_obj<LValueVector<OtherDerived>>& target) const;
        /* Getters */
        [[nodiscard]] __device__ inline ScalarType calc(size_t index) const;
        [[nodiscard]] __host__ __device__ size_t getLength() const { return mat.getDerived().getRow(); }
        [[nodiscard]] __host__ __device__ const auto& getLHS() const noexcept { return mat.getDerived(); }
        [[nodiscard]] __host__ __device__ const auto& getRHS() const noexcept { return vec.getDerived(); }
    };

    template<class MatrixType, class VectorType>
    __device__ inline typename device_obj<MatrixVectorProduct<MatrixType, VectorType>>::ScalarType
    device_obj<MatrixVectorProduct<MatrixType, VectorType>>::calc(size_t index) const {
        return mat.getDerived().row(index) * vec.getDerived();
    }

    template<class MatrixType, class VectorType>
    template<class OtherDerived>
    __device__ void device_obj<MatrixVectorProduct<MatrixType, VectorType>>::assignTo(
            device_obj<LValueVector<OtherDerived>>& target) const {
        if constexpr (MatrixOption::isColumnMatrix<MatrixType>()) {
            const auto& m = getLHS();
            const auto& v = getRHS();
            target = m.col(0).asVector() * v.calc(0);
            for (size_t i = 1; i < v.getLength(); ++i)
                target += m.col(i).asVector() * v.calc(i);
        }
        else
            Base::template assignTo<OtherDerived>(target);
    }

    template<class MatrixType, class VectorType>
    __host__ __device__ inline typename std::enable_if<MatrixType::RowAtCompile != 1, device_obj<MatrixVectorProduct<MatrixType, VectorType>>>::type
    operator*(const device_obj<RValueMatrix<MatrixType>>& mat, const device_obj<RValueVector<VectorType>>& vec) {
        assert(mat.getColumn() == vec.getLength());
        return {mat.getDerived(), vec.getDerived()};
    }
}

namespace Physica {
    template<class MatrixType, class VectorType>
    class Traits<Core::device_obj<MatrixVectorProduct<MatrixType, VectorType>>>
            : public Traits<MatrixVectorProduct<MatrixType, VectorType>> {};
}
