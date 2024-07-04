/*
 * Copyright 2024 WeiBo He.
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
    class MatrixVectorProduct<BlockMatrix<MatrixType>, VectorType>
            : public RValueVector<MatrixVectorProduct<BlockMatrix<MatrixType>, VectorType>> {
        using This = MatrixVectorProduct<BlockMatrix<MatrixType>, VectorType>;
        using Base = RValueVector<This>;
    public:
        using ScalarType = typename Base::ScalarType;
    private:
        const BlockMatrix<MatrixType>& m;
        const VectorType& v;
    public:
        MatrixVectorProduct(const BlockMatrix<MatrixType>& m_, const VectorType& v_);
        MatrixVectorProduct(const This&) = delete;
        MatrixVectorProduct(This&&) noexcept = delete;
        ~MatrixVectorProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<class OtherVector, class Executor = SequentialExecutor>
        void assignTo(LValueVector<OtherVector>& target_) const;

        [[nodiscard]] ScalarType calc(size_t index) const { throw NotImplementedException(); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        [[nodiscard]] const BlockMatrix<MatrixType>& getLHS() const noexcept { return m; }
        [[nodiscard]] const VectorType& getRHS() const noexcept { return v; }
    };

    template<class MatrixType, class VectorType>
    MatrixVectorProduct<BlockMatrix<MatrixType>, VectorType>::MatrixVectorProduct(
            const BlockMatrix<MatrixType>& m, const VectorType& v) {
        assert(m.getColumn() == v.getLength() && "[Error]: Dimensions do not match");
    }

    template<class MatrixType, class VectorType>
    template<class OtherVector, class Executor>
    void MatrixVectorProduct<BlockMatrix<MatrixType>, VectorType>::assignTo(LValueVector<OtherVector>& target_) const {
        assert(getLength() == target_.getLength() && "[Error]: Dimensions do not match");
        auto& target = target_.getDerived();
        size_t from = 0;
        for (size_t i = 0; i < m.getNumBlocks(); ++i) {
            const size_t to = m.getIndexEnds()[i];
            const auto v1 = v.segment(from, to);
            auto target1 = target.segment(from, to);
            target1 = m.getBlocks()[i] * v1;
            from = to;
        }
    }
}

namespace Physica {
    using namespace Core;

    template<class MatrixType, class VectorType>
    class Traits<MatrixVectorProduct<BlockMatrix<MatrixType>, VectorType>> {
    public:
        using ScalarType = typename MatrixType::ScalarType;
        constexpr static size_t SizeAtCompile = Dynamic;
        constexpr static size_t MaxSizeAtCompile = Dynamic;

        constexpr static bool FastAssign = true;
    };
}
