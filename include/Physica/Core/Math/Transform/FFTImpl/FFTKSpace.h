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

namespace Physica {
    template<class Derived, size_t Dim> class FFTKSpace;

    template<class Derived>
    class FFTKSpace<Derived, 1>
            : public CRTPBase<FFTKSpace<Derived, 1>>
            , public ContinuousVector<FFTKSpace<Derived, 1>> {
        using This = FFTKSpace<Derived, 1>;
        using Base = CRTPBase<This>;
        using VectorBase = ContinuousVector<This>;
    protected:
        using typename VectorBase::T;
        using typename VectorBase::Tr;
    public:
        ~FFTKSpace() = default;
        /* Operators */
        FFTKSpace& operator=(const FFTKSpace& obj);
        using VectorBase::operator=;
        using VectorBase::operator[];
        /* Operations */
        [[nodiscard]] T calc(size_t index) const;
        void invTransform(const Vector auto& data);
        [[nodiscard]] Tr parseval() const;

        using VectorBase::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        using Base::getDerived;
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Derived::rSizeToKSize(Base::getDerived().getRSpaceSize()); }
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return getLength(); }
        [[nodiscard]] __host__ __device__ size_t getRSpaceSize() const noexcept { return Base::getDerived().getRSpaceSize(); }
        [[nodiscard]] __host__ __device__ auto data(this auto&&) noexcept;
    protected:
        FFTKSpace() = default;
        FFTKSpace(const FFTKSpace&) = default;
        FFTKSpace(FFTKSpace&&) noexcept = default;
    };

    template<class Derived>
    auto FFTKSpace<Derived, 1>::operator=(const FFTKSpace<Derived, 1>& obj) -> This& {
        resize(obj.getLength());
        obj.assign(*this);
        return *this;
    }

    template<class Derived>
    auto FFTKSpace<Derived, 1>::calc(size_t index) const -> T {
        assert(index < getRSpaceSize());
        if (index < getLength())
            return (*this)[index];
        return (*this)[getRSpaceSize() - index].conjugate();
    }

    template<class Derived>
    void FFTKSpace<Derived, 1>::invTransform(const Vector auto& data) {
        VectorBase::assert_assign(data);
        *this = data;
        Base::getDerived().invTransform();
    }

    template<class Derived>
    auto FFTKSpace<Derived, 1>::parseval() const -> Tr {
        Tr result = 0;
        for (size_t i = 0; i < getRSpaceSize(); ++i)
            result += calc(i).squaredNorm();
        return result / Tr(getRSpaceSize());
    }

    template<class Derived>
    __host__ __device__ auto FFTKSpace<Derived, 1>::data(this auto&& self) noexcept {
        return self.getDerived().asComplexBuffer();
    }
    //////////////////////////////////////////////////////////////////////
    template<class Derived>
    class FFTKSpace<Derived, 2>
            : public CRTPBase<FFTKSpace<Derived, 2>>
            , public LValueMatrix<FFTKSpace<Derived, 2>> {
        using This = FFTKSpace<Derived, 2>;
        using Base = CRTPBase<This>;
        using MatrixBase = LValueMatrix<This>;
    protected:
        using typename MatrixBase::T;
    public:
        ~FFTKSpace() = default;
        /* Operators */
        FFTKSpace& operator=(const FFTKSpace& obj);
        using MatrixBase::operator=;
        using MatrixBase::operator[];
        /* Operations */
        void invTransform(const Matrix auto& data);

        using MatrixBase::resize;
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == getRow()); assert(col == getCol()); }
        /* Getters */
        using Base::getDerived;
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return Base::getDerived().getKSpaceSize()[0]; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return Base::getDerived().getKSpaceSize()[1]; }
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, size_t row, size_t col) noexcept;
    protected:
        FFTKSpace() = default;
        FFTKSpace(const FFTKSpace&) = default;
        FFTKSpace(FFTKSpace&&) noexcept = default;
    };

    template<class Derived>
    auto FFTKSpace<Derived, 2>::operator=(const FFTKSpace<Derived, 2>& obj) -> This& {
        resize(obj.getRow(), obj.getCol());
        obj.assign(*this);
        return *this;
    }

    template<class Derived>
    void FFTKSpace<Derived, 2>::invTransform(const Matrix auto& data) {
        MatrixBase::assert_assign(data);
        *this = data;
        Base::getDerived().invTransform();
    }

    template<class Derived>
    __host__ __device__ auto FFTKSpace<Derived, 2>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        assert(row < self.getRow());
        assert(col < self.getCol());
        return self.getDerived().asComplexBuffer() + row * self.getCol() + col;
    }
    //////////////////////////////////////////////////////////////////////
    template<class Derived>
    class FFTKSpace<Derived, 3>
            : public CRTPBase<FFTKSpace<Derived, 3>>
            , public LValueTensor<FFTKSpace<Derived, 3>> {
        using This = FFTKSpace<Derived, 3>;
        using Base = CRTPBase<This>;
        using TensorBase = LValueTensor<This>;
    protected:
        using typename TensorBase::T;
    public:
        ~FFTKSpace() = default;
        /* Operators */
        FFTKSpace& operator=(const FFTKSpace& obj);
        using TensorBase::operator=;
        using TensorBase::operator[];
        /* Operations */
        void invTransform(const Tensor auto& data);

        using TensorBase::resize;
        void resize([[maybe_unused]] Index3D size);
        /* Getters */
        using Base::getDerived;
        [[nodiscard]] size_t dim(int index) const noexcept { return getShape()[index]; }
        [[nodiscard]] const auto& getShape() const noexcept { return Base::getDerived().getKSpaceSize(); }
        [[nodiscard]] __host__ __device__ size_t getDimX() const noexcept { return Base::getDerived().getKSpaceSize()[0]; }
        [[nodiscard]] __host__ __device__ size_t getDimY() const noexcept { return Base::getDerived().getKSpaceSize()[1]; }
        [[nodiscard]] __host__ __device__ size_t getDimZ() const noexcept { return Base::getDerived().getKSpaceSize()[2]; }
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, Index3D index) noexcept;
    protected:
        FFTKSpace() = default;
        FFTKSpace(const FFTKSpace&) = default;
        FFTKSpace(FFTKSpace&&) noexcept = default;
    };

    template<class Derived>
    auto FFTKSpace<Derived, 3>::operator=(const FFTKSpace<Derived, 3>& obj) -> This& {
        resize(obj.getDim());
        obj.assign(*this);
        return *this;
    }

    template<class Derived>
    void FFTKSpace<Derived, 3>::invTransform(const Tensor auto& data) {
        TensorBase::assert_assign(data);
        *this = data;
        Base::getDerived().invTransform();
    }

    template<class Derived>
    void FFTKSpace<Derived, 3>::resize([[maybe_unused]] Index3D size) {
        assert(getDimX() == size[0]);
        assert(getDimY() == size[1]);
        assert(getDimZ() == size[2]);
    }

    template<class Derived>
    __host__ __device__ auto FFTKSpace<Derived, 3>::data_ptr(this auto&& self, Index3D index) noexcept {
        assert(index[0] < self.getDimX());
        assert(index[1] < self.getDimY());
        assert(index[2] < self.getDimZ());
        return self.getDerived().asComplexBuffer() + (index[0] * self.getDimY() + index[1]) * self.getDimZ() + index[2];
    }
}

namespace Physica {
    template<class T>
    class Traits<FFTKSpace<T, 1>> {
        static_assert(Traits<T>::Dim == 1, "[Error]: Inconsistent template param");
    public:
        using Derived = T;
        using ScalarType = Traits<T>::ScalarType::ComplexType;
        constexpr static size_t SizeAtCompile = Dynamic;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = true;
    };

    template<class T>
    class Traits<FFTKSpace<T, 2>> {
        static_assert(Traits<T>::Dim == 2, "[Error]: Inconsistent template param");
    public:
        using Derived = T;
        using ScalarType = Traits<T>::ScalarType::ComplexType;
        constexpr static int Major = MatrixMajor::Row;
        constexpr static size_t RowAtCompile = Dynamic;
        constexpr static size_t ColAtCompile = Dynamic;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };

    template<class T>
    class Traits<FFTKSpace<T, 3>> {
        static_assert(Traits<T>::Dim == 3, "[Error]: Inconsistent template param");
    public:
        using Derived = T;
        using ScalarType = Traits<T>::ScalarType::ComplexType;
        constexpr static int NDim = 3;
    };
}
