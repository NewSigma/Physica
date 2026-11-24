/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/StridedMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/TensorImpl/StridedTensor.h"

namespace Physica {
    template<class Derived, size_t Dim> class FFTKSpace;

    template<class Derived>
    class FFTKSpace<Derived, 1>
            : public CRTP<FFTKSpace<Derived, 1>>
            , public CompactVector<FFTKSpace<Derived, 1>> {
        using This = FFTKSpace<Derived, 1>;
        using Base = CRTP<This>;
        using VectorBase = CompactVector<This>;
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
        [[nodiscard]] __host__ __device__ auto data_handle(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Dynamic; }
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
    __host__ __device__ auto FFTKSpace<Derived, 1>::data_handle(this auto&& self) noexcept {
        return self.getDerived().asComplexBuffer();
    }
    //////////////////////////////////////////////////////////////////////
    template<class Derived>
    class FFTKSpace<Derived, 2>
            : public CRTP<FFTKSpace<Derived, 2>>
            , public StridedMatrix<FFTKSpace<Derived, 2>> {
        using This = FFTKSpace<Derived, 2>;
        using Base = CRTP<This>;
        using MatrixBase = StridedMatrix<This>;
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
        [[nodiscard]] __host__ __device__ size_t getOrder() const noexcept;
        [[nodiscard]] __host__ __device__ size_t getRowStride() const noexcept { return getCol(); }
        [[nodiscard]] __host__ __device__ auto data_handle(this auto&& self) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return MatrixMajor::Row; }
        [[nodiscard]] __host__ __device__ consteval static size_t getColStrideAtCompile() noexcept { return 1; }
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
    [[nodiscard]] __host__ __device__ size_t FFTKSpace<Derived, 2>::getOrder() const noexcept {
        assert(MatrixBase::isSquare() && "[Error]: getOrder() assumes square matrix");
        return getRow();
    }

    template<class Derived>
    __host__ __device__ auto FFTKSpace<Derived, 2>::data_handle(this auto&& self) noexcept {
        return self.getDerived().asComplexBuffer();
    }
    //////////////////////////////////////////////////////////////////////
    template<class Derived>
    class FFTKSpace<Derived, 3>
            : public CRTP<FFTKSpace<Derived, 3>>
            , public StridedTensor<FFTKSpace<Derived, 3>> {
        using This = FFTKSpace<Derived, 3>;
        using Base = CRTP<This>;
        using TensorBase = StridedTensor<This>;
    protected:
        using typename TensorBase::T;
    public:
        using typename TensorBase::ScalarType;
        using typename TensorBase::IndexType;
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
        [[nodiscard]] __host__ __device__ auto getStrides() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle(this auto&& self) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticShape() noexcept { return false; }
        [[nodiscard]] __host__ __device__ consteval static IndexType getShapeAtCompile() noexcept { return IndexType({Dynamic, Dynamic, Dynamic}); }
        [[nodiscard]] __host__ __device__ consteval static IndexType getStrideAtCompile() noexcept;
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
    __host__ __device__ auto FFTKSpace<Derived, 3>::getStrides() const noexcept {
        const size_t dimY = getDimY();
        const size_t dimZ = getDimZ();
        return IndexType({dimY * dimZ, dimZ, 1});
    }

    template<class Derived>
    __host__ __device__ auto FFTKSpace<Derived, 3>::data_handle(this auto&& self) noexcept {
        return self.getDerived().asComplexBuffer();
    }

    template<class Derived>
    __host__ __device__ consteval auto FFTKSpace<Derived, 3>::getStrideAtCompile() noexcept -> IndexType {
        return IndexType({Dynamic, Dynamic, 1});
    }
}

namespace Physica {
    template<class T>
    class Traits<FFTKSpace<T, 1>> {
        static_assert(Traits<T>::Dim == 1, "[Error]: Inconsistent template param");
    public:
        using Derived = T;
        using ScalarType = Traits<T>::ScalarType::ComplexType;
    };

    template<class T>
    class Traits<FFTKSpace<T, 2>> {
        static_assert(Traits<T>::Dim == 2, "[Error]: Inconsistent template param");
    public:
        using Derived = T;
        using ScalarType = Traits<T>::ScalarType::ComplexType;
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
