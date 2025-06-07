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

        const LazyDestroy<M> mexp;
        const LazyDestroy<V> v;
    public:
        MatExpVecProd(M&& mexp_, V&& v_);
        MatExpVecProd(const This&) = default;
        MatExpVecProd(This&&) noexcept = default;
        ~MatExpVecProd() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<Vector V1, ExecutePolicy P = Sequential>
        inline void assign(V1& target) const;

        template<Vector V1, ExecutePolicy P = Sequential>
        void assign(V1& target, Tr traceMu, ParamPair params) const;

        [[nodiscard]] ScalarType calc(size_t) const { noImpl(__func__); }
        [[nodiscard]] Tv calc_value(size_t) const { noImpl(__func__); }
        [[nodiscard]] auto calcTraceMu() const { return mexp.calcTraceMu(); }
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
    MatExpVecProd<M, V>::MatExpVecProd(M&& mexp_, V&& v_) : mexp(std::forward<M>(mexp_)), v(std::forward<V>(v_)) {
        assert(mexp.getCol() == v.getLength());
    }

    template<Matrix M, Vector V>
    template<Vector V1, ExecutePolicy P>
    inline void MatExpVecProd<M, V>::assign(V1& target) const {
        const Tr traceMu = calcTraceMu();
        assign<V1, P>(target, traceMu, calcParam<P>(traceMu));
    }

    template<Matrix M, Vector V>
    template<Vector V1, ExecutePolicy P>
    void MatExpVecProd<M, V>::assign(V1& target, Tr traceMu, ParamPair params) const {
        using BufferType = DenseVector<ScalarType, std::remove_cvref_t<V>::SizeAtCompile>;
        assert(getLength() == target.getLength());
        const Tr epsilon = std::numeric_limits<Tr>::epsilon();
        const auto& mat = mexp.getMatrix();
        const int numTaylorTerm = params.first;
        const int numSplit = params.second;
        const Tr factor = exp(traceMu / Tr(numSplit));

        BufferType buffer(getLength()), term;
        target = v;
        for (int i = 0; i < numSplit; ++i) {
            term = target;
            Tr norm1 = term.normInf();
            for (int n = 1; n <= numTaylorTerm; ++n) {
                term *= reciprocal(Tr(numSplit * n));
                {
                    auto expr = mat * term - traceMu * term;
                    expr.template assign<BufferType, P>(buffer);
                    buffer.swap(term);
                }
                const Tr norm2 = term.normInf();
                target += term;
                if (norm1 + norm2 <= epsilon * target.normInf())
                    break;
                norm1 = norm2;
            }
            target *= factor;
        }
    }

    template<Matrix M, Vector V>
    template<ExecutePolicy P>
    auto MatExpVecProd<M, V>::calcParam(Tr traceMu) const -> ParamPair {
        constexpr static Tm NormLimit = ((2 * MaxNormOrder * (MaxNormOrder + 3)) * calcTheta(MaxNumTaylorTerm)) / MaxNumTaylorTerm;
        const auto unit = UnitMatrix<ScalarType>(getLength());
        const ScalarType norm1 = (mexp.getMatrix() - unit * traceMu).template norm1_power<P>(MaxNormIteration);
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

        Trv powerNorms[MaxNormOrder];
        const Tr normalizer = reciprocal(norm1); // pow() has the risk of overflow
        for (int order = 2; order <= MaxNormOrder + 1; ++order) {
            const Tr pNorm1 = pow(mexp.getMatrix() * normalizer - (traceMu * normalizer) * unit, order).template norm1_power<P>(MaxNormIteration);
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
