/*
 * Copyright 2023-2024 Weibo He.
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
    template<class Derived, size_t Dim> class FFTKSpace;

    template<class Derived>
    class FFTKSpace<Derived, 1>
            : public CRTPBase<FFTKSpace<Derived, 1>>
            , public ContinuousVector<FFTKSpace<Derived, 1>> {
        using This = FFTKSpace<Derived, 1>;
        using Base = CRTPBase<This>;
        using VectorBase = ContinuousVector<This>;
        using RealType = typename Traits<This>::ScalarType::RealType;
    public:
        using typename VectorBase::ScalarType;
    protected:
        using typename VectorBase::PtrTy;
        using typename VectorBase::ConstPtrTy;
    public:
        ~FFTKSpace() = default;
        /* Operators */
        FFTKSpace& operator=(const FFTKSpace& obj);
        using VectorBase::operator=;
        using VectorBase::operator[];
        /* Operations */
        [[nodiscard]] ScalarType calc(size_t index) const;
        template<class VectorType>
        inline void invTransform(const RValueVector<VectorType>& data);
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Derived::rSizeToKSize(Base::getDerived().getRSpaceSize()); }
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return getLength(); }
        [[nodiscard]] __host__ __device__ size_t getRSpaceSize() const noexcept { return Base::getDerived().getRSpaceSize(); }
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t index);
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t index) const;
    protected:
        FFTKSpace() = default;
        FFTKSpace(const FFTKSpace&) = default;
        FFTKSpace(FFTKSpace&&) noexcept = default;
    };

    template<class Derived>
    FFTKSpace<Derived, 1>& FFTKSpace<Derived, 1>::operator=(const FFTKSpace<Derived, 1>& obj) {
        resize(obj.getLength());
        obj.assignTo(*this);
        return *this;
    }

    template<class Derived>
    typename FFTKSpace<Derived, 1>::ScalarType FFTKSpace<Derived, 1>::calc(size_t index) const {
        assert(index < getRSpaceSize());
        if (index < getLength())
            return (*this)[index];
        return (*this)[getRSpaceSize() - index].conjugate();
    }

    template<class Derived>
    template<class VectorType>
    inline void FFTKSpace<Derived, 1>::invTransform(const RValueVector<VectorType>& data) {
        assert(data.getLength() == getLength());
        *this = data;
        Base::getDerived().invTransform();
    }

    template<class Derived>
    __host__ __device__ inline typename FFTKSpace<Derived, 1>::PtrTy FFTKSpace<Derived, 1>::data_ptr(size_t index) {
        assert(index < getLength());
        return Base::getDerived().asComplexBuffer() + index;

    }

    template<class Derived>
    __host__ __device__ inline typename FFTKSpace<Derived, 1>::ConstPtrTy FFTKSpace<Derived, 1>::data_ptr(size_t index) const {
        return const_cast<This&>(*this).data_ptr(index);
    }
    //////////////////////////////////////////////////////////////////////
    template<class Derived>
    class FFTKSpace<Derived, 2>
            : public CRTPBase<FFTKSpace<Derived, 2>>
            , public ContinuousMatrix<FFTKSpace<Derived, 2>> {
        using This = FFTKSpace<Derived, 2>;
        using Base = CRTPBase<This>;
        using MatrixBase = ContinuousMatrix<This>;
    public:
        using typename MatrixBase::ScalarType;
    protected:
        using typename MatrixBase::PtrTy;
        using typename MatrixBase::ConstPtrTy;
    public:
        ~FFTKSpace() = default;
        /* Operators */
        FFTKSpace& operator=(const FFTKSpace& obj);
        using MatrixBase::operator=;
        using MatrixBase::operator();
        /* Operations */
        template<class MatrixType>
        inline void invTransform(const RValueMatrix<MatrixType>& data);
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == getRow()); assert(col == getCol()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return Base::getDerived().getKSpaceSize()[0]; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return Base::getDerived().getKSpaceSize()[1]; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t row, size_t col);
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t row, size_t col) const;
    protected:
        FFTKSpace() = default;
        FFTKSpace(const FFTKSpace&) = default;
        FFTKSpace(FFTKSpace&&) noexcept = default;
    };

    template<class Derived>
    FFTKSpace<Derived, 2>& FFTKSpace<Derived, 2>::operator=(const FFTKSpace<Derived, 2>& obj) {
        resize(obj.getRow(), obj.getCol());
        obj.assignTo(*this);
        return *this;
    }

    template<class Derived>
    template<class MatrixType>
    inline void FFTKSpace<Derived, 2>::invTransform(const RValueMatrix<MatrixType>& data) {
        assert(data.getRow() == getRow());
        assert(data.getCol() == getCol());
        *this = data;
        Base::getDerived().invTransform();
    }

    template<class Derived>
    inline typename FFTKSpace<Derived, 2>::PtrTy FFTKSpace<Derived, 2>::data_ptr(size_t row, size_t col) {
        assert(row < getRow());
        assert(col < getCol());
        return Base::getDerived().asComplexBuffer() + row * getCol() + col;
    }

    template<class Derived>
    inline typename FFTKSpace<Derived, 2>::ConstPtrTy FFTKSpace<Derived, 2>::data_ptr(size_t row, size_t col) const {
        return const_cast<This&>(*this).data_ptr(row, col);
    }
    //////////////////////////////////////////////////////////////////////
    template<class Derived>
    class FFTKSpace<Derived, 3>
            : public CRTPBase<FFTKSpace<Derived, 3>>
            , public LValueGrid<FFTKSpace<Derived, 3>> {
        using This = FFTKSpace<Derived, 3>;
        using Base = CRTPBase<This>;
        using GridBase = LValueGrid<This>;
    public:
        using typename GridBase::ScalarType;
        using Index3D = typename GridBase::Index3D;
    protected:
        using typename GridBase::PtrTy;
        using typename GridBase::ConstPtrTy;
    public:
        ~FFTKSpace() = default;
        /* Operators */
        FFTKSpace& operator=(const FFTKSpace& obj);
        using GridBase::operator=;
        using GridBase::operator();
        /* Operations */
        template<class GridType> inline void invTransform(const LValueGrid<GridType>& data);
        inline void resize([[maybe_unused]] Index3D size);
        using GridBase::forIndexInGrid;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getDimX() const noexcept { return Base::getDerived().getKSpaceSize()[0]; }
        [[nodiscard]] __host__ __device__ size_t getDimY() const noexcept { return Base::getDerived().getKSpaceSize()[1]; }
        [[nodiscard]] __host__ __device__ size_t getDimZ() const noexcept { return Base::getDerived().getKSpaceSize()[2]; }
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(Index3D index);
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(Index3D index) const;
    protected:
        FFTKSpace() = default;
        FFTKSpace(const FFTKSpace&) = default;
        FFTKSpace(FFTKSpace&&) noexcept = default;
    };

    template<class Derived>
    FFTKSpace<Derived, 3>& FFTKSpace<Derived, 3>::operator=(const FFTKSpace<Derived, 3>& obj) {
        resize(obj.getDim());
        obj.assignTo(*this);
        return *this;
    }

    template<class Derived>
    template<class GridType>
    inline void FFTKSpace<Derived, 3>::invTransform(const LValueGrid<GridType>& data) {
        assert(data.getDimX() == getDimX());
        assert(data.getDimY() == getDimY());
        assert(data.getDimZ() == getDimZ());
        *this = data;
        Base::getDerived().invTransform();
    }

    template<class Derived>
    inline void FFTKSpace<Derived, 3>::resize([[maybe_unused]] Index3D size) {
        assert(getDimX() == size[0]);
        assert(getDimY() == size[1]);
        assert(getDimZ() == size[2]);
    }

    template<class Derived>
    __host__ __device__ inline typename FFTKSpace<Derived, 3>::PtrTy
    FFTKSpace<Derived, 3>::data_ptr(Index3D index) {
        assert(index[0] < getDimX());
        assert(index[1] < getDimY());
        assert(index[2] < getDimZ());
        return Base::getDerived().asComplexBuffer() + (index[0] * getDimY() + index[1]) * getDimZ() + index[2];
    }

    template<class Derived>
    __host__ __device__ inline typename FFTKSpace<Derived, 3>::ConstPtrTy
    FFTKSpace<Derived, 3>::data_ptr(Index3D index) const {
        return const_cast<This&>(*this).data_ptr(index);
    }
}

namespace Physica {
    template<class T>
    class Traits<FFTKSpace<T, 1>> {
        static_assert(Traits<T>::Dim == 1, "[Error]: Inconsistent template param");
    public:
        using Derived = T;
        using ScalarType = typename Traits<T>::ScalarType::ComplexType;
        constexpr static size_t SizeAtCompile = Dynamic;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = true;
    };

    template<class T>
    class Traits<FFTKSpace<T, 2>> {
        static_assert(Traits<T>::Dim == 2, "[Error]: Inconsistent template param");
    public:
        using Derived = T;
        using ScalarType = typename Traits<T>::ScalarType::ComplexType;
        constexpr static int Option = MatrixOption::Row | MatrixOption::Element;
        constexpr static size_t RowAtCompile = Dynamic;
        constexpr static size_t ColAtCompile = Dynamic;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };

    template<class T>
    class Traits<FFTKSpace<T, 3>> {
        static_assert(Traits<T>::Dim == 3, "[Error]: Inconsistent template param");
    public:
        using Derived = T;
        using ScalarType = typename Traits<T>::ScalarType::ComplexType;
    };
}
