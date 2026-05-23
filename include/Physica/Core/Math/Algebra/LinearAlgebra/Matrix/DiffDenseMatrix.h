/*
 * Copyright 2024-2026 Weibo He.
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
    template<Scalar T, DiffMode Mode, int Order, int Major, size_t Row, size_t Col>
    class DenseMatrix<Diff<T, Mode, Order>, Major, Row, Col>
            : public CompactMatrix<DenseMatrix<Diff<T, Mode, Order> , Major, Row, Col>>
            , public std::conditional<Mode == DiffMode::Forward, CRCoro<DenseMatrix<Diff<T, Mode, Order>, Major, Row, Col>>, PlainStruct<void>>::type {
        using This = DenseMatrix<Diff<T, Mode, Order>, Major, Row, Col>;
        using Base = CompactMatrix<This>;
    public:
        using typename Base::ScalarType;
        using device_obj_type = device_obj<This>;
        using Base::isForwardDiff;
        using Base::isReverseDiff;
        template<Scalar U>
        using rebind_scalar = DenseMatrix<U, Major, Row, Col>;
    protected:
        using typename Base::Tv;
    private:
        using ValueMatrix = DenseMatrix<T, Major, Row, Col>;
        using VectorIniter = ValueMatrix::VectorIniter;
        using GradType = ScalarType::GradType;
        using GradMatrix = std::conditional<Order == 1, ValueMatrix, DenseMatrix<GradType, Major, Row, Col>>::type;

        ValueMatrix v;
        GradMatrix g;
    public:
        DenseMatrix() = default;
        explicit DenseMatrix(size_t order);
        DenseMatrix(size_t row, size_t col);
        DenseMatrix(size_t row, size_t col, T init);
        DenseMatrix(size_t row, size_t col, ScalarType init) requires(isForwardDiff());
        DenseMatrix(std::initializer_list<Tv> list);
        DenseMatrix(std::initializer_list<VectorIniter> list);
        DenseMatrix(ValueMatrix v_, GradMatrix g_);
        explicit(isReverseDiff()) DenseMatrix(const Vector auto& vec);
        explicit(isReverseDiff()) DenseMatrix(const Matrix auto& mat);
        DenseMatrix(const This&) = default;
        DenseMatrix(This&&) noexcept = default;
        ~DenseMatrix() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        using Base::operator=;
        using Base::operator[];
        /* Operations */
        void zero_grad();

        using Base::resize;
        void resize(size_t row, size_t col);
        [[nodiscard]] auto toDevice() const;
        [[nodiscard]] auto toDeviceAsync() const;
        void toDevice(device_obj<This>& obj) const;
        void toDeviceAsync(device_obj<This>& obj) const;

        template<RNG R>
        void random_uniform();
        template<RNG R>
        void random_normal();
        template<RNG R>
        void random_any(auto& distribution);

        void swap(This& __restrict obj) noexcept;
        void swap_row(size_t r1, size_t r2) noexcept;
        void swap_col(size_t c1, size_t c2) noexcept;
        /* Getters */
        [[nodiscard]] auto data(this auto&&) noexcept;
        [[nodiscard]] size_t getCol() const noexcept { return v.getCol(); }
        [[nodiscard]] size_t getRow() const noexcept { return v.getRow(); }
        [[nodiscard]] bool empty() const noexcept { return v.empty(); }

        [[nodiscard]] auto&& values(this auto&&) noexcept;
        template<int GradOrder = 1>
        [[nodiscard]] auto&& grads(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return Major; }
        [[nodiscard]] static This identity(size_t order);
        template<RNG R>
        [[nodiscard]] static auto random_uniform(size_t row, size_t col);
        template<RNG R>
        [[nodiscard]] static auto random_normal(size_t row, size_t col);
        template<RNG R>
        [[nodiscard]] static auto random_any(size_t row, size_t col, auto& distribution);
    private:
        friend class device_obj<This>;
    };
}

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order, int Major, size_t Row, size_t Col>
    class Traits<DenseMatrix<Diff<T, Mode, Order>, Major, Row, Col>> : public Traits<DenseMatrix<T, Major, Row, Col>> {
    public:
        using ScalarType = Diff<T, Mode, Order>;
    };
}

#include "MatrixImpl/DiffDenseMatrixImpl.h"
