/*
 * Copyright 2025 Weibo He.
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

#include "Physica/Core/ML/NeuralNetwork/Layer/LayerBase.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica {
    /**
     * Reference:
     * [1] ACM Transactions on Graphics, 38(5) 1-19 (2019); https://doi.org/10.1145/3341156
     */
    template<Scalar T>
    class LinearCoupler {
        using This = LinearCoupler<T>;
        using Tv = T::ValueType;
    private:
        int dim;
        int numBin;
    public:
        LinearCoupler() = default;
        LinearCoupler(int dim_, int numBin_);
        template<Scalar U>
        LinearCoupler(const LinearCoupler<U>& other);
        LinearCoupler(const This&) = default;
        LinearCoupler(This&&) noexcept = default;
        ~LinearCoupler() = default;
        /* Operators */
        This& operator=(This& obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<DNN Net>
        [[nodiscard]] CoDiff<T> forward(Net& nn, VectorND<Tv>& x) const;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] int getDim() const noexcept { return dim; }
        [[nodiscard]] int getNumBin() const noexcept { return numBin; }
    private:
        CoDiff<T> transform(const VectorND<T>& weights, VectorND<Tv>& x2) const;
        Tv calcFactor() const noexcept { return Tv(numBin) * Tv(1 - std::numeric_limits<Tv>::epsilon()); }
    };

    template<Scalar T>
    LinearCoupler<T>::LinearCoupler(int dim_, int numBin_)
            : dim(dim_), numBin(numBin_) {
        assert(dim > 0);
        assert(numBin > 0);
    }

    template<Scalar T>
    template<Scalar U>
    LinearCoupler<T>::LinearCoupler(const LinearCoupler<U>& other) : LinearCoupler(other.getDim(), other.getNumBin()) {}

    template<Scalar T>
    template<DNN Net>
    auto LinearCoupler<T>::forward(Net& nn, VectorND<Tv>& x) const -> CoDiff<T> {
        assert(getDim() == x.getLength() && "[Error]: Dimensions do not match");
        VectorND<Tv> xA(getDim());
        xA.zeros();
        xA.head(getDim() / 2) = x.head(getDim() / 2);
        VectorND<Tv> xB = x - xA;
        auto w1 = nn.forward(xA);
        auto w2 = nn.forward(xB);
        auto lnJ = transform(w1, xB) + transform(w2, xA);
        x = xA + xB;
        if constexpr (ReverseDiff<T>) {
            auto result = co_yield lnJ.value() + Tv(dim) * ln(calcFactor());
            lnJ.reverse(result.grad());
        }
        else
            co_return std::move(lnJ);
    }

    template<Scalar T>
    void LinearCoupler<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(dim, obj.dim);
        std::swap(numBin, obj.numBin);
    }

    template<Scalar T>
    auto LinearCoupler<T>::transform(const VectorND<T>& weights, VectorND<Tv>& z) const -> CoDiff<T> {
        assert(weights.getLength() == getDim() * numBin);
        const auto factor = Tv(numBin) * Tv(1 - std::numeric_limits<Tv>::epsilon());
        const auto grid = weights.reshape_col(numBin, getDim());
        const int dim = getDim();
        Array<size_t> indices(dim);
        VectorND<T> deltas(dim, 1);
        VectorND<Tv> prob(numBin);
        for (int i = 0; i < dim; ++i) {
            if (z[i].isZero())
                continue;
            const Tv tmp = z[i] * factor;
            const size_t index = tmp.toMachine();
            indices[i] = index;

            prob = softmax(grid.col(i).values());
            deltas[i] = std::max(prob[index], Tv(std::numeric_limits<Tv>::min()));

            Tv zi = tmp.mod() * prob[index];
            if (index > 0)
                zi = std::min(Tv(1), zi + prob.head(index).sum());
            z[i] = zi;
        }

        auto expr = ln(deltas);
        auto lnJ = expr.sum();
        if constexpr (ReverseDiff<T>) {
            auto tmp = co_yield lnJ.value();
            lnJ.reverse_final(tmp.grad());
            for (int i = 0; i < dim; ++i) {
                if (z[i].isZero())
                    continue;
                const auto col = grid.col(i);
                const Tv s = col.values().softmax(indices[i]);
                const Tv g = (s - square(s)) * deltas[i].grad();
                col.reverse(-g / Tv(numBin));
                col[indices[i]].reverse(g);
            }
        }
        else
            co_return std::move(lnJ);
    }
}
