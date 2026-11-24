/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/StridedMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/TensorImpl/StridedTensor.h"

namespace Physica {
    template<class Derived, size_t Dim> class FFTRSpace;

    template<class Derived>
    class FFTRSpace<Derived, 1>
            : public CRTP<FFTRSpace<Derived, 1>>
            , public CompactVector<FFTRSpace<Derived, 1>> {
        using This = FFTRSpace<Derived, 1>;
        using Base = CRTP<This>;
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
        [[nodiscard]] __host__ __device__ auto data_handle(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Dynamic; }
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
    __host__ __device__ auto FFTRSpace<Derived, 1>::data_handle(this auto&& self) noexcept {
        if constexpr (isComplex())
            return self.getDerived().asComplexBuffer();
        else
            return self.getDerived().asRealBuffer();
    }
    //////////////////////////////////////////////////////////////////////
    template<class Derived>
    class FFTRSpace<Derived, 2>
            : public CRTP<FFTRSpace<Derived, 2>>
            , public StridedMatrix<FFTRSpace<Derived, 2>> {
        using This = FFTRSpace<Derived, 2>;
        using Base = CRTP<This>;
        using MatrixBase = StridedMatrix<This>;
    public:
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
        [[nodiscard]] __host__ __device__ size_t getOrder() const noexcept;
        [[nodiscard]] __host__ __device__ size_t getRowStride() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle(this auto&& self) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return MatrixMajor::Row; }
        [[nodiscard]] __host__ __device__ consteval static size_t getColStrideAtCompile() noexcept { return 1; }
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
    __host__ __device__ size_t FFTRSpace<Derived, 2>::getOrder() const noexcept {
        assert(MatrixBase::isSquare() && "[Error]: getOrder() assumes square matrix");
        return getRow();
    }

    template<class Derived>
    __host__ __device__ size_t FFTRSpace<Derived, 2>::getRowStride() const noexcept {
        const size_t col = getCol();
        return isComplex() ? col : (col / 2 + 1) * 2;
    }

    template<class Derived>
    __host__ __device__ auto FFTRSpace<Derived, 2>::data_handle(this auto&& self) noexcept {
        if constexpr (isComplex())
            return self.getDerived().asComplexBuffer();
        else
            return self.getDerived().asRealBuffer();
    }
    //////////////////////////////////////////////////////////////////////

    template<class Derived>
    class FFTRSpace<Derived, 3>
            : public CRTP<FFTRSpace<Derived, 3>>
            , public StridedTensor<FFTRSpace<Derived, 3>> {
        using This = FFTRSpace<Derived, 3>;
        using Base = CRTP<This>;
        using TensorBase = StridedTensor<This>;
    public:
        using typename TensorBase::IndexType;
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
        [[nodiscard]] __host__ __device__ auto getStrides() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle(this auto&& self) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticShape() noexcept { return false; }
        [[nodiscard]] __host__ __device__ consteval static IndexType getShapeAtCompile() noexcept { return IndexType({Dynamic, Dynamic, Dynamic}); }
        [[nodiscard]] __host__ __device__ consteval static IndexType getStrideAtCompile() noexcept;
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
    __host__ __device__ auto FFTRSpace<Derived, 3>::getStrides() const noexcept {
        const size_t dimY = getDimY();
        const size_t dimZ = getDimZ();
        const size_t shift = isComplex() ? dimZ : (dimZ / 2 + 1) * 2;
        return IndexType({dimY * shift, shift, 1});
    }

    template<class Derived>
    __host__ __device__ auto FFTRSpace<Derived, 3>::data_handle(this auto&& self) noexcept {
        if constexpr (isComplex())
            return self.getDerived().asComplexBuffer();
        else
            return self.getDerived().asRealBuffer();
    }

    template<class Derived>
    __host__ __device__ consteval auto FFTRSpace<Derived, 3>::getStrideAtCompile() noexcept -> IndexType {
        return IndexType({Dynamic, Dynamic, 1});
    }
}

namespace Physica {
    template<class T>
    class Traits<FFTRSpace<T, 1>> {
        static_assert(Traits<T>::Dim == 1, "[Error]: Inconsistent template param");
    public:
        using Derived = T;
        using ScalarType = Traits<T>::ScalarType;
    };

    template<class T>
    class Traits<FFTRSpace<T, 2>> {
        static_assert(Traits<T>::Dim == 2, "[Error]: Inconsistent template param");
    public:
        using Derived = T;
        using ScalarType = Traits<T>::ScalarType;
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
