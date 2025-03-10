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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/UnitMatrix.h"
#include "MatrixExp.h"
#include "MatrixPow.h"

namespace Physica {
    template<Matrix T, Vector U> class MatrixVectorProduct;
    /**
     * Reference:
     * [1] SIAM J. Sci. Comput. 33(2), 488–511 (2011); https://doi.org/10.1137/100788860
     */
    template<Matrix T, Vector U>
    class MatrixVectorProduct<MatrixExp<T>, U>
            : public RValueVector<MatrixVectorProduct<MatrixExp<T>, U>> {
        using This = MatrixVectorProduct<MatrixExp<T>, U>;
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
        constexpr static bool IsFloat = ScalarType::Option == Float;
        constexpr static int MaxNumTaylorTerm = 55;
        constexpr static int MaxNormOrder = 8;
        constexpr static int MaxNormIteration = 16;
        constexpr static int BufferSize = 11;
        constexpr static float theta_single[BufferSize]{1.3E-1, 1, 2.2, 3.6, 4.9, 6.3, 7.7, 9.1, 11, 12, 13};
        constexpr static double theta_double[BufferSize]{2.4E-3, 1.4E-1, 6.4E-1, 1.4, 2.4, 3.5, 4.7, 6.0, 7.2, 8.5, 9.9};

        const MatrixExp<T>& mexp;
        const U& v;
    public:
        MatrixVectorProduct(const MatrixExp<T>& mexp_, const U& v_);
        MatrixVectorProduct(const This&) = default;
        MatrixVectorProduct(This&&) noexcept = default;
        ~MatrixVectorProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<Vector V, class Executor = SeqExecutor>
        inline void assign(V& target) const;

        template<Vector V, class Executor = SeqExecutor>
        void assign(V& target, Tr traceMu, ParamPair params) const;

        [[nodiscard]] ScalarType calc(size_t) const { noImpl(__func__); }
        [[nodiscard]] Tv calc_value(size_t) const { noImpl(__func__); }
        [[nodiscard]] auto calcTraceMu() const { return mexp.calcTraceMu(); }
        template<class Executor>
        [[nodiscard]] ParamPair calcParam(Tr traceMu) const;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        [[nodiscard]] const MatrixExp<T>& getLHS() const noexcept { return mexp; }
        [[nodiscard]] const U& getRHS() const noexcept { return v; }
    private:
        constexpr static Tm calcTheta(int numTaylorTerm);
    };

    template<Matrix T, Vector U>
    MatrixVectorProduct<MatrixExp<T>, U>::MatrixVectorProduct(
            const MatrixExp<T>& mexp_, const U& v_) : mexp(mexp_), v(v_) {
        assert(mexp.getCol() == v.getLength());
    }

    template<Matrix T, Vector U>
    template<Vector V, class Executor>
    inline void MatrixVectorProduct<MatrixExp<T>, U>::assign(V& target) const {
        const Tr traceMu = calcTraceMu();
        assign<V, Executor>(target, traceMu, calcParam<Executor>(traceMu));
    }

    template<Matrix T, Vector U>
    template<Vector V, class Executor>
    void MatrixVectorProduct<MatrixExp<T>, U>::assign(V& target, Tr traceMu, ParamPair params) const {
        using BufferType = DenseVector<ScalarType, U::SizeAtCompile>;
        assert(getLength() == target.getLength());
        const Tr epsilon = std::numeric_limits<Tr>::epsilon();
        const auto& mat = mexp.getMatrix();
        const int numTaylorTerm = params.first;
        const int numSplit = params.second;
        const Tr factor = exp(traceMu / Tr(numSplit));

        BufferType buffer, term;
        target = v;
        for (int i = 0; i < numSplit; ++i) {
            term = target;
            Tr norm1 = term.normInf();
            for (int n = 1; n <= numTaylorTerm; ++n) {
                using ExprType1 = decltype(mat * term - traceMu * term);
                term *= reciprocal(Tr(numSplit * n));
                buffer.template operator=<ExprType1, Executor>(mat * term - traceMu * term);
                buffer.swap(term);
                const Tr norm2 = term.normInf();
                target += term;
                if (norm1 + norm2 <= epsilon * target.normInf())
                    break;
                norm1 = norm2;
            }
            target *= factor;
        }
    }

    template<Matrix T, Vector U>
    template<class Executor>
    MatrixVectorProduct<MatrixExp<T>, U>::ParamPair
    MatrixVectorProduct<MatrixExp<T>, U>::calcParam(Tr traceMu) const {
        constexpr static Tm NormLimit = ((2 * MaxNormOrder * (MaxNormOrder + 3)) * calcTheta(MaxNumTaylorTerm)) / MaxNumTaylorTerm;
        const auto unit = UnitMatrix<ScalarType>(getLength());
        const ScalarType norm1 = (mexp.getMatrix() - unit * traceMu).template norm1_power<Executor>(MaxNormIteration);
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
            const Tr pNorm1 = pow(mexp.getMatrix() * normalizer - (traceMu * normalizer) * unit, order).template norm1_power<Executor>(MaxNormIteration);
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

    template<Matrix T, Vector U>
    constexpr typename MatrixVectorProduct<MatrixExp<T>, U>::Tm
    MatrixVectorProduct<MatrixExp<T>, U>::calcTheta(int numTaylorTerm) {
        assert(1 <= numTaylorTerm && numTaylorTerm <= MaxNumTaylorTerm && "[Error]: Invalid param");
        const int bufferIndex = (numTaylorTerm - 1) / 5; 
        return IsFloat ? theta_single[bufferIndex] : theta_double[bufferIndex];
    }
}

namespace Physica {
    template<Matrix T, Vector U>
    class Traits<MatrixVectorProduct<MatrixExp<T>, U>> : public Traits<MatrixVectorProduct<T, U>> {};
}
