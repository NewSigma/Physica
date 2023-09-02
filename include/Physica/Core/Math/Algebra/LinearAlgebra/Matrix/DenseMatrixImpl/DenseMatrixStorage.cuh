/*
 * Copyright 2022-2023 WeiBo He.
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

#include "DenseMatrixStorage.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.cuh"

namespace Physica::Core {
    template<class Derived>
    class device_obj<DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Element>>
            : public Utils::device_obj<Utils::Array<typename Internal::Traits<Derived>::ScalarType,
                                       Internal::Traits<Derived>::SizeAtCompile,
                                       Internal::Traits<Derived>::MaxSizeAtCompile>>
            , public Utils::CRTPBase<device_obj<Derived>, 1> {
        using T = typename Internal::Traits<Derived>::ScalarType;
        using Base = Utils::device_obj<Utils::Array<T, Internal::Traits<Derived>::SizeAtCompile, Internal::Traits<Derived>::MaxSizeAtCompile>>;
        using host_obj = DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Element>;
        using Utils::CRTPBase<device_obj<Derived>, 1>::getDerived;
    public:
        device_obj() = default;
        device_obj(size_t row, size_t column) : Base(row * column) {}
        device_obj(const host_obj& storage) : Base(storage) {}
        device_obj(const device_obj&) = default;
        device_obj(device_obj&& obj) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(device_obj other) noexcept { swap(other); return *this; }
        [[nodiscard]] __device__ T& operator()(size_t r, size_t c);
        [[nodiscard]] __device__ const T& operator()(size_t r, size_t c) const;
        /* Operations */
        [[nodiscard]] host_obj toHost() const { return host_obj(Base::toHost()); }
        void resize(size_t row, size_t column) { Base::resize(row * column); }
        using Base::swap;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return Base::getLength(); }
    };

    template<class Derived>
    __device__ typename device_obj<DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Element>>::T&
    device_obj<DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Element>>::operator()(size_t r, size_t c) {
        assert(r < getDerived().getRow() && c < getDerived().getColumn());
        return Base::operator[](getDerived().getRow() * c + r);
    }

    template<class Derived>
    __device__ const typename device_obj<DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Element>>::T&
    device_obj<DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Element>>::operator()(size_t r, size_t c) const {
        assert(r < getDerived().getRow() && c < getDerived().getColumn());
        return Base::operator[](getDerived().getRow() * c + r);
    }

    template<class Derived>
    device_obj<DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Element>>
    DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Element>::toDevice() const {
        return device_obj<This>(*this);
    }

    template<class Derived>
    class device_obj<DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Element>>
            : public Utils::device_obj<Utils::Array<typename Internal::Traits<Derived>::ScalarType,
                                       Internal::Traits<Derived>::SizeAtCompile,
                                       Internal::Traits<Derived>::MaxSizeAtCompile>>
            , public Utils::CRTPBase<device_obj<Derived>, 1> {
        using T = typename Internal::Traits<Derived>::ScalarType;
        using Base = Utils::device_obj<Utils::Array<T, Internal::Traits<Derived>::SizeAtCompile, Internal::Traits<Derived>::MaxSizeAtCompile>>;
        using host_obj = DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Element>;
        using Utils::CRTPBase<device_obj<Derived>, 1>::getDerived;
    public:
        device_obj() = default;
        device_obj(size_t row, size_t column) : Base(row * column) {}
        device_obj(const host_obj& storage) : Base(storage) {}
        device_obj(const device_obj&) = default;
        device_obj(device_obj&& obj) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(device_obj other) noexcept { swap(other); return *this; }
        [[nodiscard]] __device__ T& operator()(size_t r, size_t c);
        [[nodiscard]] __device__ const T& operator()(size_t r, size_t c) const;
        /* Operations */
        [[nodiscard]] host_obj toHost() const { return host_obj(Base::toHost()); }
        void resize(size_t row, size_t column) { Base::resize(row * column); }
        using Base::swap;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return Base::getLength(); }
    };

    template<class Derived>
    __device__ typename device_obj<DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Element>>::T&
    device_obj<DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Element>>::operator()(size_t r, size_t c) {
        assert(r < getDerived().getRow() && c < getDerived().getColumn());
        return Base::operator[](getDerived().getRow() * r + c);
    }

    template<class Derived>
    __device__ const typename device_obj<DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Element>>::T&
    device_obj<DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Element>>::operator()(size_t r, size_t c) const {
        assert(r < getDerived().getRow() && c < getDerived().getColumn());
        return Base::operator[](getDerived().getRow() * r + c);
    }

    template<class Derived>
    device_obj<DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Element>>
    DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Element>::toDevice() const {
        return device_obj<This>(*this);
    }

    template<class Derived>
    class device_obj<DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>>
            : public Utils::CRTPBase<device_obj<Derived>, 1> {
        using Base = Utils::CRTPBase<device_obj<Derived>, 1>;
        using T = typename Internal::Traits<Derived>::ScalarType;
        using host_obj = DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>;
        using VectorType = typename host_obj::VectorType;
        using ArrayType = Physica::Utils::device_obj<typename host_obj::ArrayType>;

        ArrayType array;
        size_t r;
    public:
        device_obj() : r(0) {}
        device_obj(size_t row, size_t column) : array(column, row), r(row) {}
        device_obj(const host_obj& storage) : array(storage.array), r(storage.getDerived().getRow()) {}
        device_obj(const device_obj&) = default;
        device_obj(device_obj&& obj) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(device_obj other) noexcept { swap(other); return *this; }
        [[nodiscard]] __device__ T& operator()(size_t r, size_t c);
        [[nodiscard]] __device__ const T& operator()(size_t r, size_t c) const;
        /* Operations */
        [[nodiscard]] host_obj toHost() const { return host_obj(array.toHost()); }
        void resize(size_t row, size_t column);
        void swap(device_obj& obj) noexcept { array.swap(obj); std::swap(r, obj.r); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return r * array.getLength(); }
    };

    template<class Derived>
    __device__ typename device_obj<DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>>::T&
    device_obj<DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>>::operator()(size_t r, size_t c) {
        assert(r < Base::getDerived().getRow() && c < Base::getDerived().getColumn());
        return array[c][r];
    }

    template<class Derived>
    __device__ const typename device_obj<DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>>::T&
    device_obj<DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>>::operator()(size_t r, size_t c) const {
        assert(r < Base::getDerived().getRow() && c < Base::getDerived().getColumn());
        return array[c][r];
    }

    template<class Derived>
    void device_obj<DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>>::resize(size_t row, size_t column) {
        array.resize(column);
        auto host_array = array.toHost();
        for (auto& vector : host_array)
            vector.resize(row);
        host_array.toDevice(array);
        r = row;
    }

    template<class Derived>
    device_obj<DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>>
    DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>::toDevice() const {
        return device_obj<This>(*this);
    }

    template<class Derived>
    class device_obj<DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>>
            : public Utils::CRTPBase<device_obj<Derived>, 1> {
        using Base = Utils::CRTPBase<device_obj<Derived>, 1>;
        using T = typename Internal::Traits<Derived>::ScalarType;
        using host_obj = DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>;
        using VectorType = typename host_obj::VectorType;
        using ArrayType = Physica::Utils::device_obj<typename host_obj::ArrayType>;

        ArrayType array;
        size_t c;
    public:
        device_obj() : c(0) {}
        device_obj(size_t row, size_t column) : array(row, column), c(column) {}
        device_obj(const host_obj& storage) : array(storage.array), c(storage.getDerived().getColumn()) {}
        device_obj(const device_obj&) = default;
        device_obj(device_obj&& obj) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(device_obj other) noexcept { swap(other); return *this; }
        [[nodiscard]] __device__ T& operator()(size_t r, size_t c);
        [[nodiscard]] __device__ const T& operator()(size_t r, size_t c) const;
        /* Operations */
        [[nodiscard]] host_obj toHost() const { return host_obj(array.toHost()); }
        void resize(size_t row, size_t column);
        void swap(device_obj& obj) noexcept { array.swap(obj); std::swap(c, obj.c); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return c * array.getLength(); }
    };

    template<class Derived>
    __device__ typename device_obj<DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>>::T&
    device_obj<DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>>::operator()(size_t r, size_t c) {
        assert(r < Base::getDerived().getRow() && c < Base::getDerived().getColumn());
        return array[r][c];
    }

    template<class Derived>
    __device__ const typename device_obj<DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>>::T&
    device_obj<DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>>::operator()(size_t r, size_t c) const {
        assert(r < Base::getDerived().getRow() && c < Base::getDerived().getColumn());
        return array[r][c];
    }

    template<class Derived>
    void device_obj<DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>>::resize(size_t row, size_t column) {
        array.resize(row);
        auto host_array = array.toHost();
        for (auto& vector : host_array)
            vector.resize(column);
        host_array.toDevice(array);
        c = column;
    }

    template<class Derived>
    device_obj<DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>>
    DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>::toDevice() const {
        return device_obj<This>(*this);
    }
}
