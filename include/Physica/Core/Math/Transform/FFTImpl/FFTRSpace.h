/*
 * Copyright 2023-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/TensorImpl/LValueTensor.h"

namespace Physica {
    template<class Derived, size_t Dim> class FFTRSpace;

    template<class Derived>
    class FFTRSpace<Derived, 1>
            : public CRTPBase<FFTRSpace<Derived, 1>>
            , public ContinuousVector<FFTRSpace<Derived, 1>> {
        using This = FFTRSpace<Derived, 1>;
        using Base = CRTPBase<This>;
        using VectorBase = ContinuousVector<This>;
    public:
        using typename VectorBase::ScalarType;
        using VectorBase::isComplex;
    protected:
        using typename VectorBase::PtrTy;
        using typename VectorBase::ConstPtrTy;
    public:
        ~FFTRSpace() = default;
        /* Operators */
        FFTRSpace& operator=(const FFTRSpace& obj);
        using VectorBase::operator=;
        using VectorBase::operator[];
        /* Operations */
        void transform(const Vector auto& data);
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Base::getDerived().getRSpaceSize(); }
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return getLength(); }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t index);
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t index) const;
    protected:
        FFTRSpace() = default;
        FFTRSpace(const FFTRSpace&) = default;
        FFTRSpace(FFTRSpace&&) noexcept = default;
    };

    template<class Derived>
    FFTRSpace<Derived, 1>& FFTRSpace<Derived, 1>::operator=(const FFTRSpace<Derived, 1>& obj) {
        resize(obj.getLength());
        obj.assign(*this);
        return *this;
    }

    template<class Derived>
    void FFTRSpace<Derived, 1>::transform(const Vector auto& data) {
        assert(data.getLength() == getLength());
        *this = data;
        Base::getDerived().transform();
    }

    template<class Derived>
    __host__ __device__ FFTRSpace<Derived, 1>::PtrTy FFTRSpace<Derived, 1>::data_ptr(size_t index) {
        assert(index < getLength());
        if constexpr (isComplex)
            return Base::getDerived().asComplexBuffer() + index;
        else
            return Base::getDerived().asRealBuffer() + index;
    }

    template<class Derived>
    __host__ __device__ FFTRSpace<Derived, 1>::ConstPtrTy FFTRSpace<Derived, 1>::data_ptr(size_t index) const {
        return const_cast<This&>(*this).data_ptr(index);
    }
    //////////////////////////////////////////////////////////////////////
    template<class Derived>
    class FFTRSpace<Derived, 2>
            : public CRTPBase<FFTRSpace<Derived, 2>>
            , public ContinuousMatrix<FFTRSpace<Derived, 2>> {
        using This = FFTRSpace<Derived, 2>;
        using Base = CRTPBase<This>;
        using MatrixBase = ContinuousMatrix<This>;
    public:
        using typename MatrixBase::ScalarType;
        using MatrixBase::isComplex;
    protected:
        using typename MatrixBase::PtrTy;
        using typename MatrixBase::ConstPtrTy;
    public:
        ~FFTRSpace() = default;
        /* Operators */
        FFTRSpace& operator=(const FFTRSpace& obj);
        using MatrixBase::operator=;
        using MatrixBase::operator();
        /* Operations */
        void transform(const Matrix auto& data);
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == getRow()); assert(col == getCol()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return Base::getDerived().getRSpaceSize()[0]; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return Base::getDerived().getRSpaceSize()[1]; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t row, size_t col);
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t row, size_t col) const;
    protected:
        FFTRSpace() = default;
        FFTRSpace(const FFTRSpace&) = default;
        FFTRSpace(FFTRSpace&&) noexcept = default;
    };

    template<class Derived>
    FFTRSpace<Derived, 2>& FFTRSpace<Derived, 2>::operator=(const FFTRSpace<Derived, 2>& obj) {
        resize(obj.getRow(), obj.getCol());
        obj.assign(*this);
        return *this;
    }

    template<class Derived>
    void FFTRSpace<Derived, 2>::transform(const Matrix auto& data) {
        assert(data.getRow() == getRow());
        assert(data.getCol() == getCol());
        *this = data;
        Base::getDerived().transform();
    }

    template<class Derived>
    __host__ __device__ FFTRSpace<Derived, 2>::PtrTy FFTRSpace<Derived, 2>::data_ptr(size_t row, size_t col) {
        assert(row < getRow());
        assert(col < getCol());
        const size_t numColumn = getCol();
        if constexpr (isComplex)
            return Base::getDerived().asComplexBuffer() + row * numColumn + col;
        else {
            const size_t shift = (numColumn / 2 + 1) * 2;
            return Base::getDerived().asRealBuffer() + row * shift + col;
        }
    }

    template<class Derived>
    __host__ __device__ FFTRSpace<Derived, 2>::ConstPtrTy FFTRSpace<Derived, 2>::data_ptr(size_t row, size_t col) const {
        return const_cast<This&>(*this).data_ptr(row, col);
    }
    //////////////////////////////////////////////////////////////////////
    template<class Derived>
    class FFTRSpace<Derived, 3>
            : public CRTPBase<FFTRSpace<Derived, 3>>
            , public LValueTensor<FFTRSpace<Derived, 3>> {
        using This = FFTRSpace<Derived, 3>;
        using Base = CRTPBase<This>;
        using TensorBase = LValueTensor<This>;
    public:
        using typename TensorBase::ScalarType;
        using TensorBase::isComplex;
    protected:
        using typename TensorBase::PtrTy;
        using typename TensorBase::ConstPtrTy;
    public:
        ~FFTRSpace() = default;
        /* Operators */
        FFTRSpace& operator=(const FFTRSpace& obj);
        using TensorBase::operator=;
        using TensorBase::operator();
        /* Operations */
        void transform(const Tensor auto& data);
        void resize([[maybe_unused]] Index3D size);
        /* Getters */
        [[nodiscard]] const auto& getShape() const noexcept { return Base::getDerived().getRSpaceSize(); }
        [[nodiscard]] size_t getShape(int dim) const noexcept { return getShape()[dim]; }
        [[nodiscard]] __host__ __device__ size_t getDimX() const noexcept { return Base::getDerived().getRSpaceSize()[0]; }
        [[nodiscard]] __host__ __device__ size_t getDimY() const noexcept { return Base::getDerived().getRSpaceSize()[1]; }
        [[nodiscard]] __host__ __device__ size_t getDimZ() const noexcept { return Base::getDerived().getRSpaceSize()[2]; }
        [[nodiscard]] size_t getSize() const noexcept { return TensorBase::toSize(Base::getDerived().getRSpaceSize()); }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(Index3D index);
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(Index3D index) const;
    protected:
        FFTRSpace() = default;
        FFTRSpace(const FFTRSpace&) = default;
        FFTRSpace(FFTRSpace&&) noexcept = default;
    };

    template<class Derived>
    FFTRSpace<Derived, 3>& FFTRSpace<Derived, 3>::operator=(const FFTRSpace<Derived, 3>& obj) {
        resize(obj.getDim());
        obj.assign(*this);
        return *this;
    }

    template<class Derived>
    void FFTRSpace<Derived, 3>::transform(const Tensor auto& data) {
        assert(data.getDimX() == getDimX());
        assert(data.getDimY() == getDimY());
        assert(data.getDimZ() == getDimZ());
        *this = data;
        Base::getDerived().transform();
    }

    template<class Derived>
    void FFTRSpace<Derived, 3>::resize([[maybe_unused]] Index3D size) {
        assert(getDimX() == size[0]);
        assert(getDimY() == size[1]);
        assert(getDimZ() == size[2]);
    }

    template<class Derived>
    __host__ __device__ FFTRSpace<Derived, 3>::PtrTy
    FFTRSpace<Derived, 3>::data_ptr(Index3D index) {
        assert(index[0] < getDimX());
        assert(index[1] < getDimY());
        assert(index[2] < getDimZ());
        if constexpr (isComplex)
            return Base::getDerived().asComplexBuffer() + (index[0] * getDimY() + index[1]) * getDimZ() + index[2];
        else {
            const size_t shift = (getDimZ() / 2 + 1) * 2;
            return Base::getDerived().asRealBuffer() + (index[0] * getDimY() + index[1]) * shift + index[2];
        }
    }

    template<class Derived>
    __host__ __device__ FFTRSpace<Derived, 3>::ConstPtrTy
    FFTRSpace<Derived, 3>::data_ptr(Index3D index) const {
        return const_cast<This&>(*this).data_ptr(index);
    }
}

namespace Physica {
    template<class T>
    class Traits<FFTRSpace<T, 1>> {
        static_assert(Traits<T>::Dim == 1, "[Error]: Inconsistent template param");
    public:
        using Derived = T;
        using ScalarType = Traits<T>::ScalarType;
        constexpr static size_t SizeAtCompile = Dynamic;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = true;
    };

    template<class T>
    class Traits<FFTRSpace<T, 2>> {
        static_assert(Traits<T>::Dim == 2, "[Error]: Inconsistent template param");
    public:
        using Derived = T;
        using ScalarType = Traits<T>::ScalarType;
        constexpr static int Option = MatrixOption::Row | MatrixOption::Vector;
        constexpr static size_t RowAtCompile = Dynamic;
        constexpr static size_t ColAtCompile = Dynamic;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };

    template<class T>
    class Traits<FFTRSpace<T, 3>> {
        static_assert(Traits<T>::Dim == 3, "[Error]: Inconsistent template param");
    public:
        using Derived = T;
        using ScalarType = Traits<T>::ScalarType;
        constexpr static int Dim = 3;
    };
}
