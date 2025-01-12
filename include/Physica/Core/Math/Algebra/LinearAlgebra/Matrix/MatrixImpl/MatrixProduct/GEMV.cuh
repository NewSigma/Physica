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
    template<Matrix T, Vector U>
    class device_obj<MatrixVectorProduct<T, U>>
            : public device_obj<RValueVector<MatrixVectorProduct<T, U>>> {
        static_assert(T::ColAtCompile == U::SizeAtCompile ||
                      T::ColAtCompile == Dynamic ||
                      U::SizeAtCompile == Dynamic,
                      "[Error]: Row and column do not match in matrix-vector product");
        using host_obj = MatrixVectorProduct<T, U>;
        using Base = device_obj<RValueVector<host_obj>>;
        using This = device_obj<host_obj>;
        using DeviceVector = device_obj<U>;
        using DeviceMatrix = device_obj<T>;
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
                const device_obj<T>& mat_, const device_obj<U>& vec_);
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<Vector V>
        __host__ __device__ void assign(device_obj<V>& target) const;
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

    template<Matrix T, Vector U>
    __host__ __device__ device_obj<MatrixVectorProduct<T, U>>::device_obj(
            const device_obj<T>& mat_, const device_obj<U>& vec_) {
        assert(mat_.getCol() == vec_.getLength());
        if constexpr (IsHost()) {
            mat.value = asStruct(mat_);
            vec.value = asStruct(vec_);
        }
        else {
            mat.ptr = &mat_;
            vec.ptr = &vec_;
        }
    }

    template<Matrix T, Vector U>
    template<Side Owner>
    __device__ inline device_obj<MatrixVectorProduct<T, U>>::ScalarType
    device_obj<MatrixVectorProduct<T, U>>::calc(size_t index) const {
        return getLHS<Owner>().row(index) * getRHS<Owner>();
    }

    template<Matrix T, Vector U>
    template<Vector V>
    __host__ __device__ void device_obj<MatrixVectorProduct<T, U>>::assign(device_obj<V>& target) const {
        if constexpr (IsHost())
            Base::template assign<V>(target);
        else {
            if constexpr (MatrixOption::isColMatrix<T>()) {
                const auto& m = getLHS();
                const auto& v = getRHS();
                target = m.col(0) * v.calc(0);
                for (size_t i = 1; i < v.getLength(); ++i)
                    target += m.col(i) * v.calc(i);
            }
            else
                Base::template assign<V>(target);
        }
    }

    template<Matrix T, Vector U>
    [[nodiscard]] __host__ __device__ inline std::enable_if<T::RowAtCompile != 1, device_obj<MatrixVectorProduct<T, U>>>::type
    operator*(const device_obj<T>& mat, const device_obj<U>& vec) noexcept {
        return {mat, vec};
    }
}

namespace Physica {
    template<Matrix T, Vector U>
    class Traits<Core::device_obj<MatrixVectorProduct<T, U>>> : public Traits<MatrixVectorProduct<T, U>> {};
}
