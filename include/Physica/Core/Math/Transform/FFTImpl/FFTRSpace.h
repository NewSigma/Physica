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
        template<class VectorType>
        inline void transform(const RValueVector<VectorType>& data);
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Base::getDerived().getRSpaceSize(); }
        [[nodiscard]] __host__ __device__ size_t getSize() const noexcept { return getLength(); }
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t index);
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t index) const;
    protected:
        FFTRSpace() = default;
        FFTRSpace(const FFTRSpace&) = default;
        FFTRSpace(FFTRSpace&&) noexcept = default;
    };

    template<class Derived>
    FFTRSpace<Derived, 1>& FFTRSpace<Derived, 1>::operator=(const FFTRSpace<Derived, 1>& obj) {
        resize(obj.getLength());
        obj.assignTo(*this);
        return *this;
    }

    template<class Derived>
    template<class VectorType>
    inline void FFTRSpace<Derived, 1>::transform(const RValueVector<VectorType>& data) {
        assert(data.getLength() == getLength());
        *this = data;
        Base::getDerived().transform();
    }

    template<class Derived>
    __host__ __device__ inline typename FFTRSpace<Derived, 1>::PtrTy FFTRSpace<Derived, 1>::data_ptr(size_t index) {
        assert(index < getLength());
        if constexpr (isComplex)
            return Base::getDerived().asComplexBuffer() + index;
        else
            return Base::getDerived().asRealBuffer() + index;
    }

    template<class Derived>
    __host__ __device__ inline typename FFTRSpace<Derived, 1>::ConstPtrTy FFTRSpace<Derived, 1>::data_ptr(size_t index) const {
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
        template<class MatrixType>
        inline void transform(const RValueMatrix<MatrixType>& data);
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
        obj.assignTo(*this);
        return *this;
    }

    template<class Derived>
    template<class MatrixType>
    inline void FFTRSpace<Derived, 2>::transform(const RValueMatrix<MatrixType>& data) {
        assert(data.getRow() == getRow());
        assert(data.getCol() == getCol());
        *this = data;
        Base::getDerived().transform();
    }

    template<class Derived>
    __host__ __device__ inline typename FFTRSpace<Derived, 2>::PtrTy FFTRSpace<Derived, 2>::data_ptr(size_t row, size_t col) {
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
    __host__ __device__ inline typename FFTRSpace<Derived, 2>::ConstPtrTy FFTRSpace<Derived, 2>::data_ptr(size_t row, size_t col) const {
        return const_cast<This&>(*this).data_ptr(row, col);
    }
    //////////////////////////////////////////////////////////////////////
    template<class Derived>
    class FFTRSpace<Derived, 3>
            : public CRTPBase<FFTRSpace<Derived, 3>>
            , public LValueGrid<FFTRSpace<Derived, 3>> {
        using This = FFTRSpace<Derived, 3>;
        using Base = CRTPBase<This>;
        using GridBase = LValueGrid<This>;
    public:
        using typename GridBase::ScalarType;
        using Index3D = typename GridBase::Index3D;
        using GridBase::isComplex;
    protected:
        using typename GridBase::PtrTy;
        using typename GridBase::ConstPtrTy;
    public:
        ~FFTRSpace() = default;
        /* Operators */
        FFTRSpace& operator=(const FFTRSpace& obj);
        using GridBase::operator=;
        using GridBase::operator();
        /* Operations */
        template<class GridType> inline void transform(const LValueGrid<GridType>& data);
        inline void resize([[maybe_unused]] Index3D size);
        using GridBase::forIndexInGrid;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getDimX() const noexcept { return Base::getDerived().getRSpaceSize()[0]; }
        [[nodiscard]] __host__ __device__ size_t getDimY() const noexcept { return Base::getDerived().getRSpaceSize()[1]; }
        [[nodiscard]] __host__ __device__ size_t getDimZ() const noexcept { return Base::getDerived().getRSpaceSize()[2]; }
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(Index3D index);
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(Index3D index) const;
    protected:
        FFTRSpace() = default;
        FFTRSpace(const FFTRSpace&) = default;
        FFTRSpace(FFTRSpace&&) noexcept = default;
    };

    template<class Derived>
    FFTRSpace<Derived, 3>& FFTRSpace<Derived, 3>::operator=(const FFTRSpace<Derived, 3>& obj) {
        resize(obj.getDim());
        obj.assignTo(*this);
        return *this;
    }

    template<class Derived>
    template<class GridType>
    inline void FFTRSpace<Derived, 3>::transform(const LValueGrid<GridType>& data) {
        assert(data.getDimX() == getDimX());
        assert(data.getDimY() == getDimY());
        assert(data.getDimZ() == getDimZ());
        *this = data;
        Base::getDerived().transform();
    }

    template<class Derived>
    inline void FFTRSpace<Derived, 3>::resize([[maybe_unused]] Index3D size) {
        assert(getDimX() == size[0]);
        assert(getDimY() == size[1]);
        assert(getDimZ() == size[2]);
    }

    template<class Derived>
    __host__ __device__ inline typename FFTRSpace<Derived, 3>::PtrTy
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
    __host__ __device__ inline typename FFTRSpace<Derived, 3>::ConstPtrTy
    FFTRSpace<Derived, 3>::data_ptr(Index3D index) const {
        return const_cast<This&>(*this).data_ptr(index);
    }
}

namespace Physica {
    using namespace Core;

    template<class T>
    class Traits<FFTRSpace<T, 1>> {
        static_assert(Traits<T>::Dim == 1, "[Error]: Inconsistent template param");
    public:
        using Derived = T;
        using ScalarType = typename Traits<T>::ScalarType;
        constexpr static size_t SizeAtCompile = Dynamic;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = true;
    };

    template<class T>
    class Traits<FFTRSpace<T, 2>> {
        static_assert(Traits<T>::Dim == 2, "[Error]: Inconsistent template param");
    public:
        using Derived = T;
        using ScalarType = typename Traits<T>::ScalarType;
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
        using ScalarType = typename Traits<T>::ScalarType;
    };
}
