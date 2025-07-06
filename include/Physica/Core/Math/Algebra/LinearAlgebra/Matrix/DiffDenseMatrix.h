/*
 * Copyright 2024-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.h" // IWYU pragma: export
#include "DenseMatrix.h"

namespace Physica {
    /**
     * Reference:
     * [1] Giles, M. An extended collection of matrix derivative results for forward and reverse mode algorithmic differentiation (2008); https://people.maths.ox.ac.uk/gilesm/files/NA-08-01.pdf.
     */
    template<Scalar T, DiffMode Mode, int Order, int Option, size_t Row, size_t Col>
    class DenseMatrix<Diff<T, Mode, Order>, Option, Row, Col>
            : public ContinuousMatrix<DenseMatrix<Diff<T, Mode, Order> , Option, Row, Col>>
            , public std::conditional<Mode == DiffMode::Forward, CRCoro<DenseMatrix<Diff<T, Mode, Order>, Option, Row, Col>>, PlainStruct<void>>::type {
        using This = DenseMatrix<Diff<T, Mode, Order> , Option, Row, Col>;
        using Base = ContinuousMatrix<This>;
    public:
        using typename Base::ScalarType;
        using device_obj_type = device_obj<This>;
        using Base::isForwardDiff;
        using Base::isReverseDiff;
    protected:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    private:
        using ValueMatrix = DenseMatrix<T, Option, Row, Col>;
        using initializer_list = ValueMatrix::initializer_list;
        using GradType = ScalarType::GradType;
        using GradMatrix = std::conditional<Order == 1, ValueMatrix, DenseMatrix<GradType, Option, Row, Col>>::type;

        ValueMatrix v;
        GradMatrix g;
    public:
        DenseMatrix() = default;
        DenseMatrix(size_t row, size_t col);
        DenseMatrix(size_t row, size_t col, T init);
        DenseMatrix(size_t row, size_t col, ScalarType init) requires(isForwardDiff);
        DenseMatrix(initializer_list list);
        DenseMatrix(ValueMatrix v_, GradMatrix g_);
        template<Vector V>
        explicit(isReverseDiff) DenseMatrix(const V& vec) requires(!ReverseDiff<V>);
        template<Matrix M>
        explicit(isReverseDiff) DenseMatrix(const M& mat) requires(!ReverseDiff<M>);
        DenseMatrix(const This&) = default;
        DenseMatrix(This&&) noexcept = default;
        ~DenseMatrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        using Base::operator();
        /* Operations */
        void resize(size_t row, size_t col);
        [[nodiscard]] inline auto toDevice() const;
        [[nodiscard]] inline auto toDeviceAsync() const;
        void toDevice(device_obj<This>& obj) const;
        void toDeviceAsync(device_obj<This>& obj) const;

        template<RNG R>
        inline void random_uniform();
        template<RNG R>
        inline void random_normal();
        template<RNG R>
        inline void random_any(auto& distribution);

        void swap(This& __restrict obj) noexcept;
        void swap_row(size_t r1, size_t r2) noexcept;
        void swap_col(size_t c1, size_t c2) noexcept;
        /* Getters */
        using Base::data;
        [[nodiscard]] inline PtrTy data_ptr(size_t row, size_t col) noexcept;
        [[nodiscard]] inline ConstPtrTy data_ptr(size_t row, size_t col) const noexcept;
        [[nodiscard]] size_t getCol() const noexcept { return v.getCol(); }
        [[nodiscard]] size_t getRow() const noexcept { return v.getRow(); }

        [[nodiscard]] const ValueMatrix& values() const noexcept { return v; }
        [[nodiscard]] ValueMatrix& values() noexcept { return v; }
        [[nodiscard]] const GradMatrix& grads() const noexcept { return g; }
        [[nodiscard]] GradMatrix& grads() noexcept { return g; }
        /* Static members */
        [[nodiscard]] static This unitMatrix(size_t order);
        template<RNG R>
        [[nodiscard]] inline static auto random_uniform(size_t row, size_t col);
        template<RNG R>
        [[nodiscard]] inline static auto random_normal(size_t row, size_t col);
        template<RNG R>
        [[nodiscard]] inline static auto random_any(size_t row, size_t col, auto& distribution);
    private:
        friend class device_obj<This>;
    };
}

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order, int Option, size_t Row, size_t Col>
    class Traits<DenseMatrix<Diff<T, Mode, Order>, Option, Row, Col>> : public Traits<DenseMatrix<T, Option, Row, Col>> {
    public:
        using ScalarType = Diff<T, Mode, Order>;
    };
}

#include "MatrixImpl/DiffDenseMatrixImpl.h"
