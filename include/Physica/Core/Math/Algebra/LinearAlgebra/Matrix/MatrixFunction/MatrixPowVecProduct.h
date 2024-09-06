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
    template<class MatrixType, class VectorType> class MatrixVectorProduct;

    template<class MatrixType, class VectorType>
    class MatrixVectorProduct<MatrixPow<MatrixType>, VectorType>
            : public RValueVector<MatrixVectorProduct<MatrixPow<MatrixType>, VectorType>> {
        using This = MatrixVectorProduct<MatrixPow<MatrixType>, VectorType>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
        const MatrixPow<MatrixType>& mpow;
        const VectorType& v;
    public:
        MatrixVectorProduct(const MatrixPow<MatrixType>& mpow_, const VectorType& v_);
        MatrixVectorProduct(const This&) = delete;
        MatrixVectorProduct(This&&) noexcept = delete;
        ~MatrixVectorProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<class OtherVector, class Executor = SequentialExecutor>
        inline void assignTo(LValueVector<OtherVector>& target_) const;

        [[nodiscard]] ScalarType calc(size_t) const { noImpl(); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        [[nodiscard]] const MatrixPow<MatrixType>& getLHS() const noexcept { return mpow; }
        [[nodiscard]] const VectorType& getRHS() const noexcept { return v; }
    };

    template<class MatrixType, class VectorType>
    MatrixVectorProduct<MatrixPow<MatrixType>, VectorType>::MatrixVectorProduct(
            const MatrixPow<MatrixType>& mpow_, const VectorType& v_) : mpow(mpow_), v(v_) {
        assert(mpow.getColumn() == v.getLength());
    }

    template<class MatrixType, class VectorType>
    template<class OtherVector, class Executor>
    inline void MatrixVectorProduct<MatrixPow<MatrixType>, VectorType>::assignTo(LValueVector<OtherVector>& target_) const {
        const int power = mpow.getPower();
        auto& target = target_.getDerived();
        if (power == 0) {
            target = v;
            return;
        }

        OtherVector buffer = mpow.getMatrix() * v;
        for (int i = 1; i < power; ++i) {
            buffer.swap(target);
            buffer = mpow.getMatrix() * target;
        }
        buffer.swap(target);
    }
}

namespace Physica {
    template<class MatrixType, class VectorType>
    class Traits<MatrixVectorProduct<MatrixPow<MatrixType>, VectorType>>
            : public Traits<MatrixVectorProduct<MatrixType, VectorType>> {};
}
