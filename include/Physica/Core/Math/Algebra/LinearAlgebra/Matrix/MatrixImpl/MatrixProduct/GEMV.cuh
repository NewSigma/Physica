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
        static_assert(MatrixType::ColAtCompile == VectorType::SizeAtCompile ||
                      MatrixType::ColAtCompile == Dynamic ||
                      VectorType::SizeAtCompile == Dynamic,
                      "[Error]: Row and column do not match in matrix-vector product");
        using host_obj = MatrixVectorProduct<MatrixType, VectorType>;
        using Base = device_obj<RValueVector<host_obj>>;
        using This = device_obj<host_obj>;
        using DeviceVector = device_obj<VectorType>;
        using DeviceMatrix = device_obj<MatrixType>;
    public:
        using typename Base::ScalarType;
    private:
        union {
            Physica::PlainStruct<const DeviceMatrix> value;
            const DeviceMatrix* ptr;
        } mat;

        union {
            Physica::PlainStruct<const DeviceVector> value;
            const DeviceVector* ptr;
        } vec;
    public:
        __host__ __device__ device_obj(
                const device_obj<RValueMatrix<MatrixType>>& mat_, const device_obj<RValueVector<VectorType>>& vec_);
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<class OtherDerived>
        __host__ __device__ void assignTo(device_obj<LValueVector<OtherDerived>>& target) const;
        /* Getters */
        template<Side Owner>
        [[nodiscard]] __device__ inline ScalarType calc(size_t index) const;

        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getLength() const { return getLHS<Owner>().getRow(); }

        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ const DeviceMatrix& getLHS() const noexcept {
            if constexpr (IsHost() || Owner == Side::Host)
                return mat.value.getDerived();
            else
                return *mat.ptr;
        }

        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ const DeviceVector& getRHS() const noexcept {
            if constexpr (IsHost() || Owner == Side::Host)
                return vec.value.getDerived();
            else
                return *vec.ptr;
        }
    };

    template<class MatrixType, class VectorType>
    __host__ __device__ device_obj<MatrixVectorProduct<MatrixType, VectorType>>::device_obj(
            const device_obj<RValueMatrix<MatrixType>>& mat_, const device_obj<RValueVector<VectorType>>& vec_) {
        assert(mat_.getCol() == vec_.getLength());
        if constexpr (IsHost()) {
            mat.value = asStruct(mat_.getDerived());
            vec.value = asStruct(vec_.getDerived());
        }
        else {
            mat.ptr = &mat_.getDerived();
            vec.ptr = &vec_.getDerived();
        }
    }

    template<class MatrixType, class VectorType>
    template<Side Owner>
    __device__ inline typename device_obj<MatrixVectorProduct<MatrixType, VectorType>>::ScalarType
    device_obj<MatrixVectorProduct<MatrixType, VectorType>>::calc(size_t index) const {
        return getLHS<Owner>().row(index) * getRHS<Owner>();
    }

    template<class MatrixType, class VectorType>
    template<class OtherDerived>
    __host__ __device__ void device_obj<MatrixVectorProduct<MatrixType, VectorType>>::assignTo(
            device_obj<LValueVector<OtherDerived>>& target) const {
        if constexpr (IsHost())
            Base::template assignTo<OtherDerived>(target);
        else {
            if constexpr (MatrixOption::isColMatrix<MatrixType>()) {
                const auto& m = getLHS();
                const auto& v = getRHS();
                target = m.col(0).asVector() * v.calc(0);
                for (size_t i = 1; i < v.getLength(); ++i)
                    target += m.col(i).asVector() * v.calc(i);
            }
            else
                Base::template assignTo<OtherDerived>(target);
        }
    }

    template<class MatrixType, class VectorType>
    [[nodiscard]] __host__ __device__ inline typename std::enable_if<MatrixType::RowAtCompile != 1, device_obj<MatrixVectorProduct<MatrixType, VectorType>>>::type
    operator*(const device_obj<RValueMatrix<MatrixType>>& mat, const device_obj<RValueVector<VectorType>>& vec) noexcept {
        return {mat.getDerived(), vec.getDerived()};
    }
}

namespace Physica {
    template<class MatrixType, class VectorType>
    class Traits<Core::device_obj<MatrixVectorProduct<MatrixType, VectorType>>>
            : public Traits<MatrixVectorProduct<MatrixType, VectorType>> {};
}
