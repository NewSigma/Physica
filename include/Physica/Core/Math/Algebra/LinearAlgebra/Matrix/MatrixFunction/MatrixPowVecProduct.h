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
    template<Matrix T, Vector U> class MatrixVectorProduct;

    template<Matrix T, Vector U>
    class MatrixVectorProduct<MatrixPow<T>, U>
            : public RValueVector<MatrixVectorProduct<MatrixPow<T>, U>> {
        using This = MatrixVectorProduct<MatrixPow<T>, U>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
    private:
        const MatrixPow<T>& mpow;
        const U& v;
    public:
        MatrixVectorProduct(const MatrixPow<T>& mpow_, const U& v_);
        MatrixVectorProduct(const This&) = delete;
        MatrixVectorProduct(This&&) noexcept = delete;
        ~MatrixVectorProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<LVector V, class Executor = SequentialExecutor>
        inline void assignTo(V& target) const;

        [[nodiscard]] ScalarType calc(size_t) const { noImpl("calc() is low performance and should be avoided"); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        [[nodiscard]] const MatrixPow<T>& getLHS() const noexcept { return mpow; }
        [[nodiscard]] const U& getRHS() const noexcept { return v; }
    };

    template<Matrix T, Vector U>
    MatrixVectorProduct<MatrixPow<T>, U>::MatrixVectorProduct(
            const MatrixPow<T>& mpow_, const U& v_) : mpow(mpow_), v(v_) {
        assert(mpow.getCol() == v.getLength());
    }

    template<Matrix T, Vector U>
    template<LVector V, class Executor>
    inline void MatrixVectorProduct<MatrixPow<T>, U>::assignTo(V& target) const {
        const int power = mpow.getPower();
        if (power == 0) {
            target = v;
            return;
        }

        V buffer = mpow.getMatrix() * v;
        for (int i = 1; i < power; ++i) {
            buffer.swap(target);
            buffer = mpow.getMatrix() * target;
        }
        buffer.swap(target);
    }
}

namespace Physica {
    template<Core::Matrix T, Core::Vector U>
    class Traits<Core::MatrixVectorProduct<Core::MatrixPow<T>, U>> : public Traits<Core::MatrixVectorProduct<T, U>> {};
}
