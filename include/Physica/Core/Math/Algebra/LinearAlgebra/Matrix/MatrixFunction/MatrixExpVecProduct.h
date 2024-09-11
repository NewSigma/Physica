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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/UnitMatrix.h>
#include "MatrixPow.h"

namespace Physica::Core {
    template<class MatrixType, class VectorType> class MatrixVectorProduct;
    /**
     * Reference:
     * [1] SIAM J. Sci. Comput. 33(2), 488–511 (2011); https://doi.org/10.1137/100788860
     */
    template<class MatrixType, class VectorType>
    class MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>
            : public RValueVector<MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>> {
        using This = MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>;
        using Base = RValueVector<This>;
    public:
        using typename Base::ScalarType;
        using ParamPair = std::pair<int, int>;
    private:
        using RealType = typename ScalarType::RealType;
        using TrivialType = typename RealType::TrivialType;
        constexpr static bool IsFloat = ScalarType::Option == Float;
        constexpr static int MaxNumTaylorTerm = 55;
        constexpr static int MaxNormOrder = 8;
        constexpr static int MaxNormIteration = 16;
        constexpr static int BufferSize = 11;
        constexpr static float theta_single[BufferSize]{1.3E-1, 1, 2.2, 3.6, 4.9, 6.3, 7.7, 9.1, 11, 12, 13};
        constexpr static double theta_double[BufferSize]{2.4E-3, 1.4E-1, 6.4E-1, 1.4, 2.4, 3.5, 4.7, 6.0, 7.2, 8.5, 9.9};

        const MatrixExp<MatrixType>& mexp;
        const VectorType& v;
    public:
        MatrixVectorProduct(const MatrixExp<MatrixType>& mexp_, const RValueVector<VectorType>& v_);
        MatrixVectorProduct(const This&) = delete;
        MatrixVectorProduct(This&&) noexcept = delete;
        ~MatrixVectorProduct() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<class OtherVector, class Executor = SequentialExecutor>
        inline void assignTo(LValueVector<OtherVector>& target_) const;

        template<class OtherVector, class Executor = SequentialExecutor>
        void assignTo(LValueVector<OtherVector>& target_, RealType traceMu, ParamPair params) const;

        [[nodiscard]] ScalarType calc(size_t) const { noImpl(); }
        [[nodiscard]] RealType calcTraceMu() const;
        template<class Executor>
        [[nodiscard]] ParamPair calcParam(RealType traceMu) const;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        [[nodiscard]] const MatrixExp<MatrixType>& getLHS() const noexcept { return mexp; }
        [[nodiscard]] const VectorType& getRHS() const noexcept { return v; }
    private:
        constexpr static TrivialType calcTheta(int numTaylorTerm);
    };

    template<class MatrixType, class VectorType>
    MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>::MatrixVectorProduct(
            const MatrixExp<MatrixType>& mexp_, const RValueVector<VectorType>& v_) : mexp(mexp_), v(v_.getDerived()) {
        assert(mexp.getColumn() == v.getLength());
    }

    template<class MatrixType, class VectorType>
    template<class OtherVector, class Executor>
    inline void MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>::assignTo(LValueVector<OtherVector>& target_) const {
        const RealType traceMu = calcTraceMu();
        assignTo<OtherVector, Executor>(target_, traceMu, calcParam<Executor>(traceMu));
    }

    template<class MatrixType, class VectorType>
    template<class OtherVector, class Executor>
    void MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>::assignTo(
            LValueVector<OtherVector>& target_, RealType traceMu, ParamPair params) const {
        assert(getLength() == target_.getLength());
        const RealType epsilon = std::numeric_limits<RealType>::epsilon();
        auto& target = target_.getDerived();
        const auto& mat = mexp.getMatrix();
        const int numTaylorTerm = params.first;
        const int numSplit = params.second;
        const RealType factor = exp(traceMu / RealType(numSplit));

        OtherVector buffer, term;
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

    template<class MatrixType, class VectorType>
    typename MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>::RealType
    MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>::calcTraceMu() const {
        const auto& mat = mexp.getMatrix();
        const ScalarType trace = mat.trace();
        assert(trace.getImag().isZero() && "[Error]: Not implemented");
        return trace.getReal() / RealType(getLength());
    }

    template<class MatrixType, class VectorType>
    template<class Executor>
    typename MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>::ParamPair
    MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>::calcParam(RealType traceMu) const {
        constexpr static TrivialType NormLimit = ((2 * MaxNormOrder * (MaxNormOrder + 3)) * calcTheta(MaxNumTaylorTerm)) / MaxNumTaylorTerm;
        const auto unit = UnitMatrix<ScalarType>(getLength());
        const ScalarType norm1 = (mexp.getMatrix() - unit * traceMu).template norm1_power<Executor>(MaxNormIteration);
        const bool isSmallNorm = TrivialType(norm1) <= NormLimit;
        int cost = std::numeric_limits<int>::max();
        int numMinCostTerm = 0;
        if (isSmallNorm) {
            int numSplit = 1;
            for (int numTerm = 1; numTerm <= MaxNumTaylorTerm; ++numTerm) {
                const int split = int(TrivialType(norm1) / calcTheta(numTerm)) + 1;
                const int temp = numTerm * split;
                if (cost > temp) {
                    cost = temp;
                    numMinCostTerm = numTerm;
                    numSplit = split;
                }
            }
            return std::make_pair(numMinCostTerm, numSplit);
        }

        RealType powerNorms[MaxNormOrder];
        const RealType normalizer = reciprocal(norm1); // pow() has the risk of overflow
        for (int order = 2; order <= MaxNormOrder + 1; ++order) {
            const RealType pNorm1 = pow(mexp.getMatrix() * normalizer - (traceMu * normalizer) * unit, order).template norm1_power<Executor>(MaxNormIteration);
            powerNorms[order - 2] = pow(pNorm1, reciprocal(RealType(order))) * norm1;
        }

        for (int order = 2; order <= MaxNormOrder; ++order) {
            const TrivialType powerNorm = std::max(powerNorms[order - 2], powerNorms[order - 1]).getTrivial();
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

    template<class MatrixType, class VectorType>
    constexpr typename MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>::TrivialType
    MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>::calcTheta(int numTaylorTerm) {
        assert(1 <= numTaylorTerm && numTaylorTerm <= MaxNumTaylorTerm && "[Error]: Invalid param");
        const int bufferIndex = (numTaylorTerm - 1) / 5; 
        return IsFloat ? theta_single[bufferIndex] : theta_double[bufferIndex];
    }
}

namespace Physica {
    template<class MatrixType, class VectorType>
    class Traits<Core::MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>>
            : public Traits<Core::MatrixVectorProduct<MatrixType, VectorType>> {};
}
