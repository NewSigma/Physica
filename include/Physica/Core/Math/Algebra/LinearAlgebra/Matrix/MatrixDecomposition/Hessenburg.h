/*
 * Copyright 2021-2024 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/HouseholderSequence.h"

namespace Physica::Core {
    template<Scalar T, size_t Order> class HessenburgMatrixH;
    /**
     * Decomposite matrix A like A = QHQ^H
     *
     * References:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013
     * [2] Eigen https://eigen.tuxfamily.org/
     */
    template<Scalar T, size_t Order = Dynamic>
    class Hessenburg {
        static_assert(Order > 2 || Order == Dynamic, "Unnecessary hessenburg operation on matrixes whose rank is 1 or 2");
        constexpr static size_t NormVectorLength = Order == Dynamic ? Dynamic : (Order - 2);
        using RealType = T::RealType;
        using MatrixH = HessenburgMatrixH<T, Order>;
        using HouseholderNorm = DenseVector<T, NormVectorLength>;
        using This = Hessenburg<T, Order>;
    public:
        using WorkingMatrix = DenseMatrix<T, Traits<MatrixH>::Option, Order, Order>;
    private:
        WorkingMatrix working;
        HouseholderNorm normVector;
    public:
        template<Matrix M>
        Hessenburg(const M& source);
        Hessenburg(const This&) = default;
        Hessenburg(This&&) noexcept = default;
        ~Hessenburg() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void resize(size_t size);
        void swap(This& __restrict obj);
        /* Getters */
        [[nodiscard]] size_t getSize() const noexcept { return working.getRow(); }
        [[nodiscard]] MatrixH getMatrixH() const noexcept { return MatrixH(*this); }
        [[nodiscard]] HouseholderSequence<WorkingMatrix> getMatrixQ() const noexcept;
    private:
        template<Matrix M>
        void compute(const M& source);
        friend class HessenburgMatrixH<T, Order>;
    };

    template<Scalar T, size_t Order>
    template<Matrix M>
    Hessenburg<T, Order>::Hessenburg(const M& source) {
        resize(source.getRow());
        compute(source);
    }

    template<Scalar T, size_t Order>
    void Hessenburg<T, Order>::resize(size_t size) {
        assert(size >= 3 && "[Error]: Hessenberg decomposition with a size smaller than 3 does not perform any operations");
        working.resize(size, size);
        normVector.resize(size - 2);
    }

    template<Scalar T, size_t Order>
    void Hessenburg<T, Order>::swap(This& __restrict obj) {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        working.swap(obj.working);
        normVector.swap(obj.normVector);
    }

    template<Scalar T, size_t Order>
    auto Hessenburg<T, Order>::getMatrixQ() const noexcept -> HouseholderSequence<WorkingMatrix> {
        HouseholderSequence result(working);
        result.setSize(working.getRow() - 2);
        result.setShift(1);
        return result;
    }

    template<Scalar T, size_t Order>
    template<Matrix M>
    void Hessenburg<T, Order>::compute(const M& source) {
        assert(source.getRow() == source.getCol());
        assert(source.getRow() > 2);
        const size_t order = source.getRow();
        working = source;
        for (size_t i = 0; i < order - 2; ++i) {
            auto to_col = working.col(i);
            auto temp = to_col.tail(i + 1);

            T unit;
            if (temp[0].squaredNorm() <= RealType(std::numeric_limits<RealType>::min()))
                unit = RealType(1);
            else
                unit = temp[0].unit();

            const RealType norm = householderInPlace(temp);
            normVector[i] = -norm * unit;

            if (!norm.isZero()) {
                auto target_right = working.rightCols(i + 1);
                applyHouseholder(target_right, temp);
                auto target_bottomRight = working.bottomRightCorner(i + 1);
                applyHouseholder(temp, target_bottomRight);
            }
        }
    }

    template<Scalar T, size_t Order>
    class HessenburgMatrixH : public RValueMatrix<HessenburgMatrixH<T, Order>> {
        using Base = RValueMatrix<HessenburgMatrixH<T, Order>>;
        const Hessenburg<T, Order>& hess;
    public:
        HessenburgMatrixH(const Hessenburg<T, Order>& hess_) : hess(hess_) {}
        /* Operations */
        template<Matrix M>
        void assignTo(LValueMatrix<M>& target) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return hess.getSize(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return hess.getSize(); }
    };

    template<Scalar T, size_t Order>
    template<Matrix M>
    void HessenburgMatrixH<T, Order>::assignTo(LValueMatrix<M>& target) const {
        using RealType = T::RealType;
        const size_t order = getRow();
        size_t i = 0;
        for (; i < order - 2; ++i) {
            auto fromCol = hess.working.col(i);
            auto toCol = target.col(i);
            auto copy = toCol.head(i + 1);
            copy = fromCol.head(i + 1);
            toCol[i + 1] = hess.normVector[i];
            auto zero = toCol.tail(i + 2);
            zero = RealType(0);
        }
        auto copy = target.rightCols(i);
        copy = hess.working.rightCols(i);
    }
}

namespace Physica {
    template<Core::Scalar T, size_t Order>
    class Traits<Core::HessenburgMatrixH<T, Order>> {
    public:
        using ScalarType = T;
        constexpr static int Option = Core::MatrixOption::Col | (Order == Dynamic ? Core::MatrixOption::Vector : Core::MatrixOption::Element);
        constexpr static size_t RowAtCompile = Order;
        constexpr static size_t ColAtCompile = Order;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColAtCompile;
    };
}
