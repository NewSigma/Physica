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
    namespace Internal {
        template<class T, class ReprType, class VectorType>
        class Traits<MatrixVectorProduct<Hubbard<T, ReprType>, VectorType>> {
            using MatrixType = Hubbard<T, ReprType>;
        public:
            using ScalarType = T;
            constexpr static size_t SizeAtCompile = MatrixType::RowAtCompile;
            constexpr static size_t MaxSizeAtCompile = MatrixType::MaxRowAtCompile;

            using PacketType = typename BestPacket<ScalarType, SizeAtCompile>::Type;
            constexpr static bool FastAssign = true;
        };
    }

    template<class T, class ReprType, class VectorType>
    class MatrixVectorProduct<Hubbard<T, ReprType>, VectorType>
            : public RValueVector<MatrixVectorProduct<Hubbard<T, ReprType>, VectorType>> {
        using MatrixType = Hubbard<T, ReprType>;
        using This = MatrixVectorProduct<MatrixType, VectorType>;
    public:
        using Base = RValueVector<This>;
        using typename Base::ScalarType;
        using RealType = typename ScalarType::RealType;
    private:
        const MatrixType& mat;
        const VectorType& vec;
    public:
        MatrixVectorProduct(const RValueMatrix<MatrixType>& mat_, const RValueVector<VectorType>& vec_)
                : mat(mat_.getDerived()), vec(vec_.getDerived()) {
            assert(mat.getColumn() == vec.getLength());
        }
        MatrixVectorProduct(const This&) = delete;
        MatrixVectorProduct(This&&) noexcept = delete;
        ~MatrixVectorProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<class OtherDerived>
        inline void assignTo(OtherDerived& target) const;
        /* Getters */
        [[nodiscard]] inline ScalarType calc(size_t) const { throw NotImplementedException(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const { return mat.getRow(); }
        [[nodiscard]] const MatrixType& getLHS() const noexcept { return mat; }
        [[nodiscard]] const VectorType& getRHS() const noexcept { return vec; }
    };

    template<class ScalarType, class ReprType, class VectorType>
    template<class OtherDerived>
    inline void MatrixVectorProduct<Hubbard<ScalarType, ReprType>, VectorType>::assignTo(OtherDerived& target) const {
        static_assert(std::is_base_of<LValueVector<OtherDerived>, OtherDerived>::value, "[Error]: Invalid target type");
        mat.dot(vec, target);
    }
}
