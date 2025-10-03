/*
 * Copyright 2024-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/UnitMatrix.h"
#include "MatrixExp.h"
#include "MatrixPow.h"

namespace Physica {
    /**
     * Reference:
     * [1] SIAM J. Sci. Comput. 33(2), 488–511 (2011); https://doi.org/10.1137/100788860
     */
    template<Matrix M, Vector V>
    class MatExpVecProd : public RValueVector<MatExpVecProd<M, V>> {
        using This = MatExpVecProd<M, V>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
        using ParamPair = std::pair<int, int>;
    protected:
        using typename Base::Tv;
        using typename Base::Tr;
        using typename Base::Trv;
        using Tm = Tr::MachineType;
    private:
        constexpr static bool IsFloat = ScalarType::Prec == Float;
        constexpr static int MaxNumTaylorTerm = 55;
        constexpr static int MaxNormOrder = 8;
        constexpr static int MaxNormIteration = 16;
        constexpr static int BufferSize = 11;
        constexpr static std::array<float, BufferSize> ThetaFloat32{1.3E-1, 1, 2.2, 3.6, 4.9, 6.3, 7.7, 9.1, 11, 12, 13};
        constexpr static std::array<double, BufferSize> ThetaFloat64{2.4E-3, 1.4E-1, 6.4E-1, 1.4, 2.4, 3.5, 4.7, 6.0, 7.2, 8.5, 9.9};

        LazyDestroy<M> mexp;
        LazyDestroy<V> v;
    public:
        MatExpVecProd(M mexp_, V v_);
        MatExpVecProd(const This&) = default;
        MatExpVecProd(This&&) noexcept = default;
        ~MatExpVecProd() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<bool NoFactor = false, ExecutePolicy P = Sequential>
        auto assign(Vector auto& target) const;

        template<bool NoFactor = false, ExecutePolicy P = Sequential>
        auto assign(Vector auto& target, Tr traceMu) const;

        template<bool NoFactor = false, ExecutePolicy P = Sequential>
        auto assign(Vector auto& target, Tr traceMu, ParamPair params) const -> std::conditional<NoFactor, Tr, void>::type;

        [[nodiscard]] ScalarType calc(size_t) const { noImpl(__func__); }
        [[nodiscard]] Tv calc_value(size_t) const { noImpl(__func__); }
        [[nodiscard]] Tr calcTraceMu() const { return mexp.calcTraceMu(); }
        template<ExecutePolicy P>
        [[nodiscard]] ParamPair calcParam(Tr traceMu) const;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        [[nodiscard]] const auto& getLHS() const noexcept { return mexp; }
        [[nodiscard]] const auto& getRHS() const noexcept { return v; }
    private:
        constexpr static Tm calcTheta(int numTaylorTerm);
    };

    template<Matrix M, Vector V>
    MatExpVecProd<M, V>::MatExpVecProd(M mexp_, V v_) : mexp(std::forward<M>(mexp_)), v(std::forward<V>(v_)) {
        static_assert(std::is_reference_v<M> && std::is_reference_v<V>);
        assert(mexp.getCol() == v.getLength());
    }

    template<Matrix M, Vector V>
    template<bool NoFactor, ExecutePolicy P>
    auto MatExpVecProd<M, V>::assign(Vector auto& target) const {
        return assign<NoFactor, P>(target, calcTraceMu());
    }

    template<Matrix M, Vector V>
    template<bool NoFactor, ExecutePolicy P>
    auto MatExpVecProd<M, V>::assign(Vector auto& target, Tr traceMu) const {
        return assign<NoFactor, P>(target, traceMu, calcParam<P>(traceMu));
    }

    template<Matrix M, Vector V>
    template<bool NoFactor, ExecutePolicy P>
    auto MatExpVecProd<M, V>::assign(Vector auto& target, Tr traceMu, ParamPair params) const -> std::conditional<NoFactor, Tr, void>::type {
        constexpr size_t SizeAtCompile1 = std::remove_cvref_t<decltype(target)>::SizeAtCompile;
        constexpr size_t SizeAtCompile2 = std::max(SizeAtCompile1, Base::SizeAtCompile);
        using BufferType = DenseVector<ScalarType, SizeAtCompile2>;

        assert(getLength() == target.getLength());
        const Tr epsilon = std::numeric_limits<Tr>::epsilon();
        const auto& mat = mexp.getMatrix();
        const int numTaylorTerm = params.first;
        const int numSplit = params.second;
        const Tr factor = exp(traceMu / Tr(numSplit));

        BufferType buffer(getLength());
        BufferType term(getLength());
        Tr lnNorm1 = 0;
        v.assign(target);
        for (int i = 0; i < numSplit; ++i) {
            Tr norm1 = target.normInf();
            if constexpr (NoFactor) {
                if (norm1.isZero()) [[unlikely]] {
                    target.zeros();
                    break;
                }
                target *= reciprocal(norm1); // Avoid overflow
                lnNorm1 += ln(norm1);
            }

            target.assign(term);
            for (int n = 1; n <= numTaylorTerm; ++n) {
                term *= reciprocal(Tr(numSplit * n));
                {
                    auto expr = mat * term - traceMu * term;
                    expr.template assign<P>(buffer);
                    buffer.swap(term);
                }
                const Tr norm2 = term.normInf();
                target += term;
                if (norm1 + norm2 <= epsilon * target.normInf())
                    break;
                norm1 = norm2;
            }

            if constexpr (!NoFactor)
                target *= factor;
        }

        if constexpr (NoFactor)
            return lnNorm1;
    }

    template<Matrix M, Vector V>
    template<ExecutePolicy P>
    auto MatExpVecProd<M, V>::calcParam(Tr traceMu) const -> ParamPair {
        constexpr static Tm NormLimit = ((2 * MaxNormOrder * (MaxNormOrder + 3)) * calcTheta(MaxNumTaylorTerm)) / MaxNumTaylorTerm;
        const auto matI = UnitMatrix<Trv>(getLength());
        const Tr norm1 = (mexp.getMatrix() - matI * traceMu).template norm1_power<P>(MaxNormIteration);
        const bool isSmallNorm = Tm(norm1) <= NormLimit;
        int cost = std::numeric_limits<int>::max();
        int numMinCostTerm = 0;
        if (isSmallNorm) {
            int numSplit = 1;
            for (int numTerm = 1; numTerm <= MaxNumTaylorTerm; ++numTerm) {
                const int split = int(Tm(norm1) / calcTheta(numTerm)) + 1;
                const int temp = numTerm * split;
                if (cost > temp) {
                    cost = temp;
                    numMinCostTerm = numTerm;
                    numSplit = split;
                }
            }
            return std::make_pair(numMinCostTerm, numSplit);
        }

        Array<Trv, MaxNormOrder> powerNorms{};
        const Tr normalizer = reciprocal(norm1); // pow() has the risk of overflow
        for (int order = 2; order <= MaxNormOrder + 1; ++order) {
            const Tr pNorm1 = pow(mexp.getMatrix() * normalizer - (traceMu * normalizer) * matI, order).template norm1_power<P>(MaxNormIteration);
            powerNorms[order - 2] = pow(pNorm1.value(), reciprocal(Trv(order))) * norm1.value();
        }

        for (int order = 2; order <= MaxNormOrder; ++order) {
            const Tm powerNorm = std::max(powerNorms[order - 2], powerNorms[order - 1]).toMachine();
            for (int numTerm = order * (order - 1) - 1; numTerm <= MaxNumTaylorTerm; ++numTerm) {
                const int temp = numTerm * int(powerNorm / calcTheta(numTerm)) + numTerm;
                if (cost > temp) {
                    cost = temp;
                    numMinCostTerm = numTerm;
                }
            }
        }
        return std::make_pair(numMinCostTerm, std::max(cost / numMinCostTerm, 1));
    }

    template<Matrix M, Vector V>
    constexpr auto MatExpVecProd<M, V>::calcTheta(int numTaylorTerm) -> Tm {
        assert(1 <= numTaylorTerm && numTaylorTerm <= MaxNumTaylorTerm && "[Error]: Invalid param");
        const int bufferIndex = (numTaylorTerm - 1) / 5; 
        return IsFloat ? ThetaFloat32[bufferIndex] : ThetaFloat64[bufferIndex];
    }
}

namespace Physica {
    template<Matrix M, Vector V>
    class Traits<MatExpVecProd<M, V>> : public Traits<GEMV<M, V>> {};
}
