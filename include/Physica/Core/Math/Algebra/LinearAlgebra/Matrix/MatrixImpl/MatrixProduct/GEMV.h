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

namespace Physica::Core {
    template<class MatrixType, class VectorType>
    class MatrixVectorProduct : public RValueVector<MatrixVectorProduct<MatrixType, VectorType>> {
        using This = MatrixVectorProduct<MatrixType, VectorType>;
    public:
        using Base = RValueVector<This>;
        using typename Base::ScalarType;
    private:
        const MatrixType& mat;
        const VectorType& vec;
    public:
        MatrixVectorProduct(const RValueMatrix<MatrixType>& mat_, const RValueVector<VectorType>& vec_)
                : mat(mat_.getDerived()), vec(vec_.getDerived()) {
            assert(mat.getCol() == vec.getLength());
        }
        MatrixVectorProduct(const This&) = delete;
        MatrixVectorProduct(This&&) noexcept = delete;
        ~MatrixVectorProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<class OtherDerived, class Executor = SequentialExecutor>
        inline void assignTo(LValueVector<OtherDerived>& target) const;
        /* Getters */
        [[nodiscard]] inline ScalarType calc(size_t index) const;
        [[nodiscard]] __host__ __device__ size_t getLength() const { return mat.getRow(); }
        [[nodiscard]] const MatrixType& getLHS() const noexcept { return mat; }
        [[nodiscard]] const VectorType& getRHS() const noexcept { return vec; }
    };

    template<class MatrixType, class VectorType>
    template<class OtherDerived, class Executor>
    inline void MatrixVectorProduct<MatrixType, VectorType>::assignTo(LValueVector<OtherDerived>& target) const {
        if constexpr (MatrixOption::isColMatrix<MatrixType>()) {
            target = mat.col(0).asVector() * vec.calc(0);
            for (size_t i = 1; i < vec.getLength(); ++i)
                target += mat.col(i).asVector() * vec.calc(i);
        }
        else {
            for (size_t i = 0; i < getLength(); ++i)
                target[i] = calc(i);
            
            constexpr bool isContinuous = std::is_base_of<ContinuousVector<OtherDerived>, OtherDerived>::value;
            if constexpr (isContinuous && Base::isReverseDiff)
                target.getDerived().makeContinuous();
        }
    }

    template<class MatrixType, class VectorType>
    inline typename MatrixVectorProduct<MatrixType, VectorType>::ScalarType MatrixVectorProduct<MatrixType, VectorType>::calc(size_t index) const {
        return mat.row(index) * vec;
    }

    template<class MatrixType, class VectorType>
    [[nodiscard]] inline typename std::enable_if<MatrixType::RowAtCompile != 1, MatrixVectorProduct<MatrixType, VectorType>>::type
    operator*(const RValueMatrix<MatrixType>& mat, const RValueVector<VectorType>& vec) noexcept {
        assert(mat.getCol() == vec.getLength());
        return MatrixVectorProduct(mat.getDerived(), vec.getDerived());
    }

    template<class MatrixType, class VectorType>
    [[nodiscard]] inline typename std::enable_if<MatrixType::RowAtCompile == 1 && MatrixType::ColAtCompile == 1,
                                                 typename Internal::BinaryScalarOpReturnType<typename MatrixType::ScalarType,
                                                                                             typename VectorType::ScalarType>::Type>::type
    operator*(const RValueMatrix<MatrixType>& mat, const RValueVector<VectorType>& vec) {
        assert(mat.getCol() == vec.getLength());
        return mat.row(0) * vec;
    }
}

namespace Physica {
    template<class MatrixType, class VectorType>
    class Traits<Core::MatrixVectorProduct<MatrixType, VectorType>> {
        static_assert(MatrixType::ColAtCompile == VectorType::SizeAtCompile ||
                      MatrixType::ColAtCompile == Dynamic ||
                      VectorType::SizeAtCompile == Dynamic,
                      "Row and column do not match in matrix product");
    public:
        using ScalarType = typename Core::Internal::BinaryScalarOpReturnType<typename MatrixType::ScalarType,
                                                                             typename VectorType::ScalarType>::Type;
        constexpr static size_t SizeAtCompile = MatrixType::RowAtCompile;
        constexpr static bool FastAssign = Core::MatrixOption::isColMatrix<MatrixType>();
        constexpr static bool FastPacket = false;
    };
}
