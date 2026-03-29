/*
 * Copyright 2022-2026 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseHermiteMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/HouseholderSequence.h"

namespace Physica {
    template<Scalar T, size_t Order> class TridiagonalMatrixT;
    /**
     * Decomposite hermite matrix A like A = QTQ^H
     * 
     * References:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013:458-459
     */
    template<Scalar T, size_t Order = Dynamic>
    class Tridiagonalization {
        static_assert(Order > 2 || Order == Dynamic,
                      "Unnecessary hessenburg operation on matrixes whose rank is 1 or 2");
        using This = Tridiagonalization<T, Order>;
        using RealType = T::RealType;
        using MatrixT = TridiagonalMatrixT<T, Order>;
        using SymmMatrix = DenseSymmMatrix<T, Order>;
        using HermiteMatrix = DenseHermiteMatrix<T, Order>;

        constexpr static size_t NormVectorLength = Order == Dynamic ? Dynamic : (Order - 2);
        constexpr static size_t BufferLength = Order == Dynamic ? Dynamic : (Order - 1);
        using HouseholderNorm = DenseVector<T, NormVectorLength>;
        using BufferVector = DenseVector<T, BufferLength>;
    public:
        using WorkingMatrix = std::conditional<T::isComplex, HermiteMatrix, SymmMatrix>::type;
    private:
        WorkingMatrix working;
        HouseholderNorm normBuffer;
        BufferVector buffer;
    public:
        Tridiagonalization(size_t size);
        Tridiagonalization(const Matrix auto& source);
        Tridiagonalization(const This&) = default;
        Tridiagonalization(This&&) noexcept = default;
        ~Tridiagonalization() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void compute(const Matrix auto& source);
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] MatrixT getMatrixT() const noexcept { return MatrixT(*this); }
        [[nodiscard]] HouseholderSequence<WorkingMatrix> getMatrixQ() const noexcept;

        friend class TridiagonalMatrixT<T, Order>;
    };

    template<Scalar T, size_t Order>
    Tridiagonalization<T, Order>::Tridiagonalization(size_t size)
            : working(size, size)
            , normBuffer(size - 2)
            , buffer(size - 1) {
        assert(size > 2 && "[Error]: Matrix is too small, no need for tridiagonalization");
    }

    template<Scalar T, size_t Order>
    Tridiagonalization<T, Order>::Tridiagonalization(const Matrix auto& source)
            : Tridiagonalization(source.getRow()) {
        compute(source);
    }

    template<Scalar T, size_t Order>
    void Tridiagonalization<T, Order>::compute(const Matrix auto& source) {
        const size_t order = source.getRow();
        working = source;
        for (size_t i = 0; i < order - 2; ++i) {
            auto to_col = working.col(i);
            auto temp = to_col.tail(i + 1);

            T factor;
            if (temp.calc(0).squaredNorm() <= RealType(std::numeric_limits<RealType>::min()))
                factor = T(1);
            else
                factor = unit(temp.calc_value(0));

            auto p = buffer.head(order - i - 1);
            const RealType norm = temp.householder(p);
            for (size_t j = 0; j < p.getLength(); ++j)
                working[i, i + j + 1] = p[j].conjugate();
            normBuffer[i] = -norm * factor;

            if (!norm.isZero()) {
                const T factor = working[i, i + 1];
                working[i, i + 1] = T(1);
                auto corner = working.bottomRightCorner(i + 1);
                p = factor * (corner * temp);
                p -= (p.conjugate() * temp * factor * T(0.5)) * temp;
                for (size_t r = 0; r < corner.getRow(); ++r)
                    for (size_t c = r; c < corner.getCol(); ++c)
                        working[r + i + 1, c + i + 1] -= temp.calc(r) * p[c].conjugate() + temp.calc(c).conjugate() * T(p[r]);
                working[i, i + 1] = factor;
            }
        }
    }

    template<Scalar T, size_t Order>
    void Tridiagonalization<T, Order>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        working.swap(obj.working);
        normBuffer.swap(obj.normBuffer);
        buffer.swap(obj.buffer);
    }

    template<Scalar T, size_t Order>
    HouseholderSequence<typename Tridiagonalization<T, Order>::WorkingMatrix>
    Tridiagonalization<T, Order>::getMatrixQ() const noexcept {
        HouseholderSequence result(working);
        result.setSize(working.getRow() - 2);
        result.setShift(1);
        return result;
    }

    template<Scalar T, size_t Order>
    class TridiagonalMatrixT : public RValueMatrix<TridiagonalMatrixT<T, Order>> {
        using This = TridiagonalMatrixT<T, Order>;
        using RealType = T::RealType;

        const Tridiagonalization<T, Order>& tri;
    public:
        TridiagonalMatrixT(const Tridiagonalization<T, Order>& tri_);
        TridiagonalMatrixT(const This&) = delete;
        TridiagonalMatrixT(This&&) noexcept = delete;
        ~TridiagonalMatrixT() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto& target) const;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return tri.working.getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return tri.working.getCol(); }
    };

    template<Scalar T, size_t Order>
    TridiagonalMatrixT<T, Order>::TridiagonalMatrixT(const Tridiagonalization<T, Order>& tri_) : tri(tri_) {}

    template<Scalar T, size_t Order>
    void TridiagonalMatrixT<T, Order>::assign(Matrix auto& target) const {
        const size_t order = getRow();
        target = RealType(0);
        target[0, 0] = tri.working.calc(0, 0);
        target[1, 0] = tri.normBuffer[0];
        size_t i = 1;
        for (; i < order - 2; ++i) {
            target[i - 1, i] = target[i, i - 1].conjugate();
            target[i, i] = tri.working.calc(i, i);
            target[i + 1, i] = tri.normBuffer[i];
        }
        target[i - 1, i] = target[i, i - 1].conjugate();
        target[i, i] = tri.working.calc(i, i);
        target[i + 1, i] = tri.working.calc(i + 1, i);
        ++i;
        target[i - 1, i] = target[i, i - 1].conjugate();
        target[i, i] = tri.working.calc(i, i);
    }
}

namespace Physica {
    template<Scalar T, size_t Order>
    class Traits<TridiagonalMatrixT<T, Order>> {
    public:
        using ScalarType = T;
        constexpr static int Major = MatrixMajor::BothMajor;
        constexpr static size_t RowAtCompile = Order;
        constexpr static size_t ColAtCompile = Order;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}
