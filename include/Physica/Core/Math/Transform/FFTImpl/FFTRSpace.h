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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/TensorImpl/LValueTensor.h"

namespace Physica {
    template<class Derived, size_t Dim> class FFTRSpace;

    template<class Derived>
    class FFTRSpace<Derived, 1>
            : public CRTPBase<FFTRSpace<Derived, 1>>
            , public CompactVector<FFTRSpace<Derived, 1>> {
        using This = FFTRSpace<Derived, 1>;
        using Base = CRTPBase<This>;
        using VectorBase = CompactVector<This>;
    public:
        using typename VectorBase::ScalarType;
        using VectorBase::isComplex;
    public:
        ~FFTRSpace() = default;
        /* Operators */
        FFTRSpace& operator=(const FFTRSpace& obj);
        using VectorBase::operator=;
        using VectorBase::operator[];
        /* Operations */
        void transform(const Vector auto& data);

        using VectorBase::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        using Base::getDerived;
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Base::getDerived().getRSpaceSize(); }
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return getLength(); }
        [[nodiscard]] __host__ __device__ auto data(this auto&&) noexcept;
    protected:
        FFTRSpace() = default;
        FFTRSpace(const FFTRSpace&) = default;
        FFTRSpace(FFTRSpace&&) noexcept = default;
    };

    template<class Derived>
    auto FFTRSpace<Derived, 1>::operator=(const FFTRSpace<Derived, 1>& obj) -> This& {
        resize(obj.getLength());
        obj.assign(*this);
        return *this;
    }

    template<class Derived>
    void FFTRSpace<Derived, 1>::transform(const Vector auto& data) {
        VectorBase::assert_assign(data);
        *this = data;
        Base::getDerived().transform();
    }

    template<class Derived>
    __host__ __device__ auto FFTRSpace<Derived, 1>::data(this auto&& self) noexcept {
        if constexpr (isComplex())
            return self.getDerived().asComplexBuffer();
        else
            return self.getDerived().asRealBuffer();
    }
    //////////////////////////////////////////////////////////////////////
    template<class Derived>
    class FFTRSpace<Derived, 2>
            : public CRTPBase<FFTRSpace<Derived, 2>>
            , public LValueMatrix<FFTRSpace<Derived, 2>> {
        using This = FFTRSpace<Derived, 2>;
        using Base = CRTPBase<This>;
        using MatrixBase = LValueMatrix<This>;
    public:
        using typename MatrixBase::ScalarType;
        using MatrixBase::isComplex;
    public:
        ~FFTRSpace() = default;
        /* Operators */
        FFTRSpace& operator=(const FFTRSpace& obj);
        using MatrixBase::operator=;
        using MatrixBase::operator[];
        /* Operations */
        void transform(const Matrix auto& data);

        using MatrixBase::resize;
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == getRow()); assert(col == getCol()); }
        /* Getters */
        using Base::getDerived;
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return Base::getDerived().getRSpaceSize()[0]; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return Base::getDerived().getRSpaceSize()[1]; }
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&& self, size_t row, size_t col) noexcept;
    protected:
        FFTRSpace() = default;
        FFTRSpace(const FFTRSpace&) = default;
        FFTRSpace(FFTRSpace&&) noexcept = default;
    };

    template<class Derived>
    auto FFTRSpace<Derived, 2>::operator=(const FFTRSpace<Derived, 2>& obj) -> This& {
        resize(obj.getRow(), obj.getCol());
        obj.assign(*this);
        return *this;
    }

    template<class Derived>
    void FFTRSpace<Derived, 2>::transform(const Matrix auto& data) {
        MatrixBase::assert_assign(data);
        *this = data;
        Base::getDerived().transform();
    }

    template<class Derived>
    __host__ __device__ auto FFTRSpace<Derived, 2>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        assert(row < self.getRow());
        assert(col < self.getCol());
        const size_t numColumn = self.getCol();
        if constexpr (isComplex())
            return self.getDerived().asComplexBuffer() + row * numColumn + col;
        else {
            const size_t shift = (numColumn / 2 + 1) * 2;
            return self.getDerived().asRealBuffer() + row * shift + col;
        }
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
    public:
        ~FFTRSpace() = default;
        /* Operators */
        FFTRSpace& operator=(const FFTRSpace& obj);
        using TensorBase::operator=;
        using TensorBase::operator[];
        /* Operations */
        void transform(const Tensor auto& data);

        using TensorBase::resize;
        void resize([[maybe_unused]] Index3D size);
        /* Getters */
        using Base::getDerived;
        [[nodiscard]] size_t dim(int index) const noexcept { return getShape()[index]; }
        [[nodiscard]] const auto& getShape() const noexcept { return Base::getDerived().getRSpaceSize(); }
        [[nodiscard]] __host__ __device__ size_t getDimX() const noexcept { return Base::getDerived().getRSpaceSize()[0]; }
        [[nodiscard]] __host__ __device__ size_t getDimY() const noexcept { return Base::getDerived().getRSpaceSize()[1]; }
        [[nodiscard]] __host__ __device__ size_t getDimZ() const noexcept { return Base::getDerived().getRSpaceSize()[2]; }
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, Index3D index) noexcept;
    protected:
        FFTRSpace() = default;
        FFTRSpace(const FFTRSpace&) = default;
        FFTRSpace(FFTRSpace&&) noexcept = default;
    };

    template<class Derived>
    auto FFTRSpace<Derived, 3>::operator=(const FFTRSpace<Derived, 3>& obj) -> This& {
        resize(obj.getDim());
        obj.assign(*this);
        return *this;
    }

    template<class Derived>
    void FFTRSpace<Derived, 3>::transform(const Tensor auto& data) {
        TensorBase::assert_assign(data);
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
    __host__ __device__ auto FFTRSpace<Derived, 3>::data_ptr(this auto&& self, Index3D index) noexcept {
        assert(index[0] < self.getDimX());
        assert(index[1] < self.getDimY());
        assert(index[2] < self.getDimZ());
        if constexpr (isComplex())
            return self.getDerived().asComplexBuffer() + (index[0] * self.getDimY() + index[1]) * self.getDimZ() + index[2];
        else {
            const size_t shift = (self.getDimZ() / 2 + 1) * 2;
            return self.getDerived().asRealBuffer() + (index[0] * self.getDimY() + index[1]) * shift + index[2];
        }
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
    };

    template<class T>
    class Traits<FFTRSpace<T, 2>> {
        static_assert(Traits<T>::Dim == 2, "[Error]: Inconsistent template param");
    public:
        using Derived = T;
        using ScalarType = Traits<T>::ScalarType;
        constexpr static int Major = MatrixMajor::Row;
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
        constexpr static int NDim = 3;
    };
}
