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

#include <Physica/Core/MultiPrecision/Diff.h>
#include "DenseMatrix.h"

namespace Physica::Core {
    template<class T, int Order, int Option>
    class DenseMatrix<Diff<T, DiffMode::Forward, Order>, Option> : public ContinuousMatrix<DenseMatrix<Diff<T, DiffMode::Forward, Order> , Option>> {
        using This = DenseMatrix<Diff<T, DiffMode::Forward, Order> , Option>;
        using Base = ContinuousMatrix<This>;
    public:
        using typename Base::ScalarType;
    protected:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    private:
        using ValueMatrix = DenseMatrix<T, Option>;
        using GradType = typename ScalarType::GradType;
        using GradMatrix = typename std::conditional<Order == 1, ValueMatrix, DenseMatrix<GradType, Option>>::type;

        ValueMatrix values;
        GradMatrix grads;
    public:
        DenseMatrix() = default;
        DenseMatrix(size_t row, size_t col);
        DenseMatrix(size_t row, size_t col, T init);
        template<class OtherMatrix>
        DenseMatrix(const RValueMatrix<OtherMatrix>& mat);
        template<class VectorType>
        DenseMatrix(const RValueVector<VectorType>& vec);
        DenseMatrix(const This&) = default;
        DenseMatrix(This&&) noexcept = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        using Base::operator();
        /* Operations */
        void resize(size_t row, size_t col);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        using Base::data;
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t row, size_t col) noexcept;
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t row, size_t col) const noexcept;
        [[nodiscard]] size_t getColumn() const noexcept { return values.getColumn(); }
        [[nodiscard]] size_t getRow() const noexcept { return values.getRow(); }
    };

    template<class PlainScalar, int Option, int Order>
    class Diff<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>
            : public RValueMatrix<Diff<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>>
            , public DenseMatrixDim<Diff<DenseMatrix<PlainScalar, Option>, DiffMode::Reverse, Order>, Dynamic, Dynamic> {
        using PlainMatrix = DenseMatrix<PlainScalar, Option>;
        using This = Diff<PlainMatrix, DiffMode::Reverse, Order>;
        using Base = RValueMatrix<This>;
        using Dim = DenseMatrixDim<This, Dynamic, Dynamic>;
        using TracerType = DiffTracer<PlainScalar, Order>;
        using SegmentType = typename TracerType::SegmentType;
    public:
        using device_obj_type = device_obj<This>;
        using ScalarType = typename Base::ScalarType;
    private:
        SegmentType& traceSeg;
    public:
        Diff();
        Diff(size_t row, size_t column);
        Diff(PlainMatrix values);
        Diff(const Diff&) = default;
        Diff(Diff&&) noexcept = default;
        ~Diff() = default;
        /* Operators */
        Diff& operator=(const Diff&) = default;
        Diff& operator=(Diff&&) noexcept = default;
        /* Operations */
        [[nodiscard]] inline ScalarType calc(size_t row, size_t col) const;

        template<class RandomGenerator>
        inline void random_uniform(RandomGenerator& gen);
        template<class RandomGenerator>
        inline void random_normal(RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        inline void random_any(Distribution& dist, RandomGenerator& gen);
        void swap(Diff& obj) noexcept { std::swap(*this, obj); }
        /* Getters */
        using Dim::getColumn;
        using Dim::getRow;
        [[nodiscard]] size_t getSize() const noexcept { return traceSeg.getLength(); }
        /* Static members */
        template<class RandomGenerator>
        [[nodiscard]] inline static This random_uniform(size_t row, size_t column, RandomGenerator& gen);
        template<class RandomGenerator>
        [[nodiscard]] inline static This random_normal(size_t row, size_t column, RandomGenerator& gen);
        template<class Distribution, class RandomGenerator>
        [[nodiscard]] inline static This random_any(size_t row, size_t column, Distribution& dist, RandomGenerator& gen);
    private:
        friend class device_obj<This>;
    };
}

namespace Physica {
    template<class T, int Order, int Option>
    class Traits<Core::DenseMatrix<Core::Diff<T, Core::DiffMode::Forward, Order>, Option>> : public Traits<Core::DenseMatrix<T, Option>> {
    public:
        using ScalarType = Core::Diff<T, Core::DiffMode::Forward, Order>;
    };

    template<class T, int Option, int Order>
    class Traits<Core::Diff<Core::DenseMatrix<T, Option>, Core::DiffMode::Reverse, Order>> : public Traits<Core::DenseMatrix<T, Option>> {
        static_assert(!T::isDifferentiable, "[Error]: Nested Diff<> is not allowed");
    public:
        using ScalarType = Core::Diff<T, Core::DiffMode::Reverse, Order>;
    };
}

#include "DiffDenseMatrixImpl/DiffDenseMatrixImpl.h"
