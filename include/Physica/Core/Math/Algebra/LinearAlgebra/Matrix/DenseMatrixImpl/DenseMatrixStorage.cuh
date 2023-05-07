/*
 * Copyright 2022 WeiBo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.cuh"

namespace Physica::Core::Internal {
    template<class Derived>
    class device_obj<DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Element>>
            : public device_obj<Utils::Array<typename Traits<Derived>::ScalarType,
                                             Traits<Derived>::SizeAtCompile,
                                             Traits<Derived>::MaxSizeAtCompile>>
            , public Utils::CRTPBase<device_obj<Derived>, 1> {
        using T = typename Traits<Derived>::ScalarType;
        using Base = device_obj<Utils::Array<T, Traits<Derived>::SizeAtCompile, Traits<Derived>::MaxSizeAtCompile>>;
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
            : public device_obj<Utils::Array<typename Traits<Derived>::ScalarType,
                                             Traits<Derived>::SizeAtCompile,
                                             Traits<Derived>::MaxSizeAtCompile>>
            , public Utils::CRTPBase<device_obj<Derived>, 1> {
        using T = typename Traits<Derived>::ScalarType;
        using Base = device_obj<Utils::Array<T, Traits<Derived>::SizeAtCompile, Traits<Derived>::MaxSizeAtCompile>>;
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
            : public Physica::Utils::device_obj<Utils::Array<Vector<typename Traits<Derived>::ScalarType, Traits<Derived>::RowAtCompile, Traits<Derived>::MaxRowAtCompile>,
                                                             Traits<Derived>::ColumnAtCompile,
                                                             Traits<Derived>::MaxColumnAtCompile>>
            , public Utils::CRTPBase<device_obj<Derived>, 1> {
        using T = typename Traits<Derived>::ScalarType;
        using VectorType = Vector<T, Traits<Derived>::RowAtCompile, Traits<Derived>::MaxRowAtCompile>;
        using Base = Physica::Utils::device_obj<Utils::Array<VectorType, Traits<Derived>::ColumnAtCompile, Traits<Derived>::MaxColumnAtCompile>>;
        using host_obj = DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>;
        using Utils::CRTPBase<device_obj<Derived>, 1>::getDerived;

        size_t r;
    public:
        device_obj() : r(0) {}
        device_obj(size_t row, size_t column) : Base(column, row), r(row) {}
        device_obj(const host_obj& storage) : Base(storage), r(storage.getDerived().getRow()) {}
        device_obj(const device_obj&) = default;
        device_obj(device_obj&& obj) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(device_obj other) noexcept { swap(other); return *this; }
        [[nodiscard]] __device__ T& operator()(size_t r, size_t c);
        [[nodiscard]] __device__ const T& operator()(size_t r, size_t c) const;
        /* Operations */
        [[nodiscard]] host_obj toHost() const { return host_obj(Base::toHost()); }
        void swap(device_obj& obj) noexcept { Base::swap(obj); std::swap(r, obj.r); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return r * Base::getLength(); }
    };

    template<class Derived>
    __device__ typename device_obj<DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>>::T&
    device_obj<DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>>::operator()(size_t r, size_t c) {
        assert(r < getDerived().getRow() && c < getDerived().getColumn());
        return Base::operator[](c)[r];
    }

    template<class Derived>
    __device__ const typename device_obj<DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>>::T&
    device_obj<DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>>::operator()(size_t r, size_t c) const {
        assert(r < getDerived().getRow() && c < getDerived().getColumn());
        return Base::operator[](c)[r];
    }

    template<class Derived>
    device_obj<DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>>
    DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>::toDevice() const {
        return device_obj<This>(*this);
    }

    template<class Derived>
    class device_obj<DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>>
            : public device_obj<Utils::Array<Vector<typename Traits<Derived>::ScalarType, Traits<Derived>::ColumnAtCompile, Traits<Derived>::MaxColumnAtCompile>,
                                                    Traits<Derived>::RowAtCompile,
                                                    Traits<Derived>::MaxRowAtCompile>>
            , public Utils::CRTPBase<device_obj<Derived>, 1> {
        using T = typename Traits<Derived>::ScalarType;
        using VectorType = Vector<T, Traits<Derived>::ColumnAtCompile, Traits<Derived>::MaxColumnAtCompile>;
        using Base = device_obj<Utils::Array<VectorType, Traits<Derived>::RowAtCompile, Traits<Derived>::MaxRowAtCompile>>;
        using host_obj = DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>;
        using Utils::CRTPBase<device_obj<Derived>, 1>::getDerived;

        size_t c;
    public:
        device_obj() : c(0) {}
        device_obj(size_t row, size_t column) : Base(row, column), c(column) {}
        device_obj(const host_obj& storage) : Base(storage), c(storage.getDerived().getColumn()) {}
        device_obj(const device_obj&) = default;
        device_obj(device_obj&& obj) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(device_obj other) noexcept { swap(other); return *this; }
        [[nodiscard]] __device__ T& operator()(size_t r, size_t c);
        [[nodiscard]] __device__ const T& operator()(size_t r, size_t c) const;
        /* Operations */
        [[nodiscard]] host_obj toHost() const { return host_obj(Base::toHost()); }
        void swap(device_obj& obj) noexcept { Base::swap(obj); std::swap(c, obj.c); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return c * Base::getLength(); }
    };

    template<class Derived>
    __device__ typename device_obj<DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>>::T&
    device_obj<DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>>::operator()(size_t r, size_t c) {
        assert(r < getDerived().getRow() && c < getDerived().getColumn());
        return Base::operator[](r)[c];
    }

    template<class Derived>
    __device__ const typename device_obj<DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>>::T&
    device_obj<DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>>::operator()(size_t r, size_t c) const {
        assert(r < getDerived().getRow() && c < getDerived().getColumn());
        return Base::operator[](r)[c];
    }

    template<class Derived>
    device_obj<DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>>
    DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>::toDevice() const {
        return device_obj<This>(*this);
    }
}
