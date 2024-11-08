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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/HouseholderSequence.h>

namespace Physica::Core {
    template<class ScalarType, size_t Order> class HessenburgMatrixH;
    /**
     * Decomposite matrix A like A = QHQ^H
     *
     * References:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013
     * [2] Eigen https://eigen.tuxfamily.org/
     */
    template<class ScalarType, size_t Order = Dynamic>
    class Hessenburg {
        static_assert(Order > 2 || Order == Dynamic, "Unnecessary hessenburg operation on matrixes whose rank is 1 or 2");
        constexpr static size_t NormVectorLength = Order == Dynamic ? Dynamic : (Order - 2);
        using RealType = typename ScalarType::RealType;
        using MatrixH = HessenburgMatrixH<ScalarType, Order>;
        using HouseholderNorm = Vector<ScalarType, NormVectorLength>;
        using This = Hessenburg<ScalarType, Order>;
    public:
        using WorkingMatrix = DenseMatrix<ScalarType, Traits<MatrixH>::Option, Order, Order>;
    private:
        WorkingMatrix working;
        HouseholderNorm normVector;
    public:
        template<class MatrixType>
        Hessenburg(const RValueMatrix<MatrixType>& source);
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
        template<class MatrixType>
        void compute(const MatrixType& source);
        friend class HessenburgMatrixH<ScalarType, Order>;
    };

    template<class ScalarType, size_t Order>
    template<class MatrixType>
    Hessenburg<ScalarType, Order>::Hessenburg(const RValueMatrix<MatrixType>& source) {
        resize(source.getRow());
        compute(source.getDerived());
    }

    template<class ScalarType, size_t Order>
    void Hessenburg<ScalarType, Order>::resize(size_t size) {
        assert(size >= 3 && "[Error]: Hessenberg decomposition with a size smaller than 3 does not perform any operations");
        working.resize(size, size);
        normVector.resize(size - 2);
    }

    template<class ScalarType, size_t Order>
    void Hessenburg<ScalarType, Order>::swap(This& __restrict obj) {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        working.swap(obj.working);
        normVector.swap(obj.normVector);
    }

    template<class ScalarType, size_t Order>
    HouseholderSequence<typename Hessenburg<ScalarType, Order>::WorkingMatrix>
    Hessenburg<ScalarType, Order>::getMatrixQ() const noexcept {
        HouseholderSequence result(working);
        result.setSize(working.getRow() - 2);
        result.setShift(1);
        return result;
    }

    template<class ScalarType, size_t Order>
    template<class MatrixType>
    void Hessenburg<ScalarType, Order>::compute(const MatrixType& source) {
        assert(source.getRow() == source.getCol());
        assert(source.getRow() > 2);
        const size_t order = source.getRow();
        working = source;
        for (size_t i = 0; i < order - 2; ++i) {
            auto to_col = working.col(i);
            auto temp = to_col.tail(i + 1);

            ScalarType unit;
            if (temp[0].squaredNorm() <= RealType(std::numeric_limits<RealType>::min()))
                unit = ScalarType(1);
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

    template<class ScalarType, size_t Order>
    class HessenburgMatrixH : public RValueMatrix<HessenburgMatrixH<ScalarType, Order>> {
        using Base = RValueMatrix<HessenburgMatrixH<ScalarType, Order>>;
        const Hessenburg<ScalarType, Order>& hess;
    public:
        HessenburgMatrixH(const Hessenburg<ScalarType, Order>& hess_) : hess(hess_) {}
        /* Operations */
        template<class OtherMatrix>
        void assignTo(LValueMatrix<OtherMatrix>& target) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return hess.getSize(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return hess.getSize(); }
    };

    template<class ScalarType, size_t Order>
    template<class OtherMatrix>
    void HessenburgMatrixH<ScalarType, Order>::assignTo(LValueMatrix<OtherMatrix>& target) const {
        const size_t order = getRow();
        size_t i = 0;
        for (; i < order - 2; ++i) {
            auto fromCol = hess.working.col(i);
            auto toCol = target.col(i);
            auto copy = toCol.head(i + 1);
            copy = fromCol.head(i + 1);
            toCol[i + 1] = hess.normVector[i];
            auto zero = toCol.tail(i + 2);
            zero = ScalarType(0);
        }
        auto copy = target.rightCols(i);
        copy = hess.working.rightCols(i);
    }
}

namespace Physica {
    using namespace Core;

    template<class T, size_t Order>
    class Traits<HessenburgMatrixH<T, Order>> {
    public:
        using ScalarType = T;
        constexpr static int Option = MatrixOption::Col | (Order == Dynamic ? MatrixOption::Vector : MatrixOption::Element);
        constexpr static size_t RowAtCompile = Order;
        constexpr static size_t ColumnAtCompile = Order;
        constexpr static size_t SizeAtCompile = RowAtCompile * ColumnAtCompile;
    };
}
