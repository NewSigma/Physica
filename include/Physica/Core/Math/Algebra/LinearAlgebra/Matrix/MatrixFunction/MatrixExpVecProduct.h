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
#include "MatrixPow.h"

namespace Physica::Core {
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
    private:
        using RealType = ScalarType::RealType;
        using RealValue = RealType::ValueType;
        using MachineType = RealType::MachineType;
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
        MatrixVectorProduct(const This&) = delete;
        MatrixVectorProduct(This&&) noexcept = delete;
        ~MatrixVectorProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<LVector V, class Executor = SequentialExecutor>
        inline void assignTo(V& target) const;

        template<LVector V, class Executor = SequentialExecutor>
        void assignTo(V& target, RealType traceMu, ParamPair params) const;

        [[nodiscard]] ScalarType calc(size_t) const { noImpl("calc() is low performance and should be avoided"); }
        [[nodiscard]] RealType calcTraceMu() const;
        template<class Executor>
        [[nodiscard]] ParamPair calcParam(RealType traceMu) const;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        [[nodiscard]] const MatrixExp<T>& getLHS() const noexcept { return mexp; }
        [[nodiscard]] const U& getRHS() const noexcept { return v; }
    private:
        constexpr static MachineType calcTheta(int numTaylorTerm);
    };

    template<Matrix T, Vector U>
    MatrixVectorProduct<MatrixExp<T>, U>::MatrixVectorProduct(
            const MatrixExp<T>& mexp_, const U& v_) : mexp(mexp_), v(v_) {
        assert(mexp.getCol() == v.getLength());
    }

    template<Matrix T, Vector U>
    template<LVector V, class Executor>
    inline void MatrixVectorProduct<MatrixExp<T>, U>::assignTo(V& target) const {
        const RealType traceMu = calcTraceMu();
        assignTo<V, Executor>(target, traceMu, calcParam<Executor>(traceMu));
    }

    template<Matrix T, Vector U>
    template<LVector V, class Executor>
    void MatrixVectorProduct<MatrixExp<T>, U>::assignTo(V& target, RealType traceMu, ParamPair params) const {
        using BufferType = DenseVector<ScalarType, U::SizeAtCompile>;
        assert(getLength() == target.getLength());
        const RealType epsilon = std::numeric_limits<RealType>::epsilon();
        const auto& mat = mexp.getMatrix();
        const int numTaylorTerm = params.first;
        const int numSplit = params.second;
        const RealType factor = exp(traceMu / RealType(numSplit));

        BufferType buffer, term;
        target = v;
        for (int i = 0; i < numSplit; ++i) {
            term = target;
            RealType norm1 = term.normInf();
            for (int n = 1; n <= numTaylorTerm; ++n) {
                using ExprType1 = decltype(mat * term - traceMu * term);
                term *= reciprocal(RealType(numSplit * n));
                buffer.template operator=<ExprType1, Executor>(mat * term - traceMu * term);
                buffer.swap(term);
                const RealType norm2 = term.normInf();
                target += term;
                if (norm1 + norm2 <= epsilon * target.normInf())
                    break;
                norm1 = norm2;
            }
            target *= factor;
        }
    }

    template<Matrix T, Vector U>
    MatrixVectorProduct<MatrixExp<T>, U>::RealType
    MatrixVectorProduct<MatrixExp<T>, U>::calcTraceMu() const {
        const auto& mat = mexp.getMatrix();
        const ScalarType trace = mat.trace();
        assert(trace.imag().isZero() && "[Error]: Not implemented");
        return trace.real() / RealType(getLength());
    }

    template<Matrix T, Vector U>
    template<class Executor>
    MatrixVectorProduct<MatrixExp<T>, U>::ParamPair
    MatrixVectorProduct<MatrixExp<T>, U>::calcParam(RealType traceMu) const {
        constexpr static MachineType NormLimit = ((2 * MaxNormOrder * (MaxNormOrder + 3)) * calcTheta(MaxNumTaylorTerm)) / MaxNumTaylorTerm;
        const auto unit = UnitMatrix<ScalarType>(getLength());
        const ScalarType norm1 = (mexp.getMatrix() - unit * traceMu).template norm1_power<Executor>(MaxNormIteration);
        const bool isSmallNorm = MachineType(norm1) <= NormLimit;
        int cost = std::numeric_limits<int>::max();
        int numMinCostTerm = 0;
        if (isSmallNorm) {
            int numSplit = 1;
            for (int numTerm = 1; numTerm <= MaxNumTaylorTerm; ++numTerm) {
                const int split = int(MachineType(norm1) / calcTheta(numTerm)) + 1;
                const int temp = numTerm * split;
                if (cost > temp) {
                    cost = temp;
                    numMinCostTerm = numTerm;
                    numSplit = split;
                }
            }
            return std::make_pair(numMinCostTerm, numSplit);
        }

        RealValue powerNorms[MaxNormOrder];
        const RealType normalizer = reciprocal(norm1); // pow() has the risk of overflow
        for (int order = 2; order <= MaxNormOrder + 1; ++order) {
            const RealType pNorm1 = pow(mexp.getMatrix() * normalizer - (traceMu * normalizer) * unit, order).template norm1_power<Executor>(MaxNormIteration);
            powerNorms[order - 2] = pow(pNorm1.getValue(), reciprocal(RealValue(order))) * norm1.getValue();
        }

        for (int order = 2; order <= MaxNormOrder; ++order) {
            const MachineType powerNorm = std::max(powerNorms[order - 2], powerNorms[order - 1]).toMachine();
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
    constexpr typename MatrixVectorProduct<MatrixExp<T>, U>::MachineType
    MatrixVectorProduct<MatrixExp<T>, U>::calcTheta(int numTaylorTerm) {
        assert(1 <= numTaylorTerm && numTaylorTerm <= MaxNumTaylorTerm && "[Error]: Invalid param");
        const int bufferIndex = (numTaylorTerm - 1) / 5; 
        return IsFloat ? theta_single[bufferIndex] : theta_double[bufferIndex];
    }
}

namespace Physica {
    template<Matrix T, Vector U>
    class Traits<Core::MatrixVectorProduct<MatrixExp<T>, U>>
            : public Traits<Core::MatrixVectorProduct<T, U>> {};
}
