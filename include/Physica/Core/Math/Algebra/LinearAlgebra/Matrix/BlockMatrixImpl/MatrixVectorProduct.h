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
    template<Matrix T, Vector U>
    class MatrixVectorProduct<BlockMatrix<T>, U>
            : public RValueVector<MatrixVectorProduct<BlockMatrix<T>, U>> {
        using This = MatrixVectorProduct<BlockMatrix<T>, U>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    private:
        const BlockMatrix<T>& m;
        const U& v;
    public:
        MatrixVectorProduct(const BlockMatrix<T>& m_, const U& v_);
        MatrixVectorProduct(const This&) = delete;
        MatrixVectorProduct(This&&) noexcept = delete;
        ~MatrixVectorProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<Vector V, class Executor = SequentialExecutor>
        void assignTo(V& target_) const;

        [[nodiscard]] ScalarType calc(size_t) const { noImpl(__func__); }
        [[nodiscard]] ValueType calc_value(size_t) const { noImpl(__func__); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        [[nodiscard]] const BlockMatrix<T>& getLHS() const noexcept { return m; }
        [[nodiscard]] const U& getRHS() const noexcept { return v; }
    };

    template<Matrix T, Vector U>
    MatrixVectorProduct<BlockMatrix<T>, U>::MatrixVectorProduct(
            const BlockMatrix<T>& m, const U& v) {
        assert(m.getCol() == v.getLength() && "[Error]: Dimensions do not match");
    }

    template<Matrix T, Vector U>
    template<Vector V, class Executor>
    void MatrixVectorProduct<BlockMatrix<T>, U>::assignTo(V& target) const {
        assert(getLength() == target.getLength() && "[Error]: Dimensions do not match");
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
    template<Matrix T, Vector U>
    class Traits<MatrixVectorProduct<BlockMatrix<T>, U>> {
    public:
        using ScalarType = T::ScalarType;
        constexpr static size_t SizeAtCompile = Dynamic;

        constexpr static bool FastAssign = true;
        constexpr static bool FastPacket = false;
    };
}
