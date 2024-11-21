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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.h"
#include "DenseMatrix.h"

namespace Physica::Core {
    template<Scalar T, int Order, int Option, size_t Row, size_t Col>
    class DenseMatrix<Diff<T, DiffMode::Forward, Order>, Option, Row, Col> : public ContinuousMatrix<DenseMatrix<Diff<T, DiffMode::Forward, Order> , Option, Row, Col>> {
        using This = DenseMatrix<Diff<T, DiffMode::Forward, Order> , Option, Row, Col>;
        using Base = ContinuousMatrix<This>;
    public:
        using typename Base::ScalarType;
    protected:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    private:
        using ValueMatrix = DenseMatrix<T, Option, Row, Col>;
        using GradType = ScalarType::GradType;
        using GradMatrix = std::conditional<Order == 1, ValueMatrix, DenseMatrix<GradType, Option, Row, Col>>::type;

        ValueMatrix values;
        GradMatrix grads;
    public:
        DenseMatrix() = default;
        DenseMatrix(size_t row, size_t col);
        DenseMatrix(size_t row, size_t col, T init);
        template<Matrix M>
        DenseMatrix(const M& mat);
        template<Vector V>
        DenseMatrix(const V& vec);
        DenseMatrix(const This&) = default;
        DenseMatrix(This&&) noexcept = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        using Base::operator();
        /* Operations */
        void resize(size_t row, size_t col);
        void swap(This& __restrict obj) noexcept;
        void swap_row(size_t r1, size_t r2) noexcept;
        void swap_col(size_t c1, size_t c2) noexcept;
        /* Getters */
        using Base::data;
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t row, size_t col) noexcept;
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t row, size_t col) const noexcept;
        [[nodiscard]] size_t getCol() const noexcept { return values.getCol(); }
        [[nodiscard]] size_t getRow() const noexcept { return values.getRow(); }
    };

    template<Scalar T, int Option, int Order>
    class Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>
            : public RValueMatrix<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>>
            , public DenseMatrixDim<Diff<DenseMatrix<T, Option>, DiffMode::Reverse, Order>, Dynamic, Dynamic> {
        using PlainMatrix = DenseMatrix<T, Option>;
        using This = Diff<PlainMatrix, DiffMode::Reverse, Order>;
        using Base = RValueMatrix<This>;
        using Dim = DenseMatrixDim<This, Dynamic, Dynamic>;
        using TracerType = DiffTracer<T, Order>;
        using SegmentType = TracerType::SegmentType;
    public:
        using device_obj_type = device_obj<This>;
        using ScalarType = Base::ScalarType;
    private:
        SegmentType& traceSeg;
    public:
        Diff();
        Diff(size_t row, size_t col);
        Diff(PlainMatrix values);
        Diff(const Diff&) = default;
        Diff(Diff&&) noexcept = default;
        ~Diff() = default;
        /* Operators */
        Diff& operator=(const Diff&) = default;
        Diff& operator=(Diff&&) noexcept = default;
        /* Operations */
        [[nodiscard]] inline ScalarType calc(size_t row, size_t col) const;

        template<class RandomType>
        inline void random_uniform(RandomType& gen);
        template<class RandomType>
        inline void random_normal(RandomType& gen);
        template<class Distribution, class RandomType>
        inline void random_any(Distribution& dist, RandomType& gen);
        void swap(Diff& obj) noexcept { std::swap(*this, obj); }
        /* Getters */
        using Dim::getCol;
        using Dim::getRow;
        [[nodiscard]] size_t getSize() const noexcept { return traceSeg.getLength(); }
        /* Static members */
        template<class RandomType>
        [[nodiscard]] inline static This random_uniform(size_t row, size_t col, RandomType& gen);
        template<class RandomType>
        [[nodiscard]] inline static This random_normal(size_t row, size_t col, RandomType& gen);
        template<class Distribution, class RandomType>
        [[nodiscard]] inline static This random_any(size_t row, size_t col, Distribution& dist, RandomType& gen);
    private:
        friend class device_obj<This>;
    };
}

namespace Physica {
    template<Scalar T, int Order, int Option, size_t Row, size_t Col>
    class Traits<Core::DenseMatrix<Core::Diff<T, Core::DiffMode::Forward, Order>, Option, Row, Col>> : public Traits<Core::DenseMatrix<T, Option>> {
    public:
        using ScalarType = Core::Diff<T, Core::DiffMode::Forward, Order>;
    };

    template<Scalar T, int Option, int Order>
    class Traits<Core::Diff<Core::DenseMatrix<T, Option>, Core::DiffMode::Reverse, Order>> : public Traits<Core::DenseMatrix<T, Option>> {
        static_assert(!T::isDifferentiable, "[Error]: Nested Diff<> is not allowed");
    public:
        using ScalarType = Core::Diff<T, Core::DiffMode::Reverse, Order>;
    };
}

#include "DiffDenseMatrixImpl/DiffDenseMatrixImpl.h"
