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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/UnitMatrix.h>
#include "MatrixPow.h"

namespace Physica::Core {
    template<class MatrixType, class VectorType> class MatrixVectorProduct;

    namespace Internal {
        template<class MatrixType, class VectorType>
        class Traits<MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>>
                : public Traits<MatrixVectorProduct<MatrixType, VectorType>> {};
    }
    /**
     * Reference:
     * [1] SIAM J. Sci. Comput. 33(2), 488–511 (2011); https://doi.org/10.1137/100788860
     */
    template<class MatrixType, class VectorType>
    class MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>
            : public RValueVector<MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>> {
        using This = MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>;
        using Base = RValueVector<This>;

        constexpr static int BufferSize = 11;
        constexpr static float theta_single[BufferSize]{1.3E-1, 1, 2.2, 3.6, 4.9, 6.3, 7.7, 9.1, 11, 12, 13};
        constexpr static double theta_double[BufferSize]{2.4E-3, 1.4E-1, 6.4E-1, 1.4, 2.4, 3.5, 4.7, 6.0, 7.2, 8.5, 9.9};
    public:
        using typename Base::ScalarType;
    private:
        using RealType = typename ScalarType::RealType;
        using TrivialType = typename RealType::TrivialType;
        constexpr static bool IsFloat = ScalarType::Option == Float;
        constexpr static int MaxNumTaylorTerm = 55;
        constexpr static int MaxNormOrder = 8;
        constexpr static int MaxNormIteration = 16;

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
        void assignTo(LValueVector<OtherVector>& target_) const;

        [[nodiscard]] ScalarType calc(size_t) const { throw NotImplementedException(); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return v.getLength(); }
        [[nodiscard]] const MatrixExp<MatrixType>& getLHS() const noexcept { return mexp; }
        [[nodiscard]] const VectorType& getRHS() const noexcept { return v; }
    private:
        constexpr static TrivialType calcTheta(int numTaylorTerm);
        std::pair<int, int> calcParam(RealType traceMu) const;
        inline RealType calcPowerNorm(RealType traceMu, int order) const;
    };

    template<class MatrixType, class VectorType>
    MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>::MatrixVectorProduct(
            const MatrixExp<MatrixType>& mexp_, const RValueVector<VectorType>& v_) : mexp(mexp_), v(v_.getDerived()) {
        assert(mexp.getColumn() == v.getLength());
    }

    template<class MatrixType, class VectorType>
    template<class OtherVector, class Executor>
    void MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>::assignTo(LValueVector<OtherVector>& target_) const {
        assert(getLength() == target_.getLength());
        const RealType epsilon = std::numeric_limits<RealType>::epsilon();
        const size_t length = getLength();
        auto& target = target_.getDerived();

        const auto& mat = mexp.getMatrix();
        const ScalarType trace = mat.trace();
        assert(trace.getImag().isZero() && "[Error]: Not implemented");
        const RealType traceMu = trace.getReal() / RealType(length);

        const auto pair = calcParam(traceMu);
        const int numTaylorTerm = pair.first;
        const int numSplit = pair.second;
        const RealType factor = exp(trace.getReal() / RealType(length * numSplit));

        OtherVector buffer, term;
        target = v;
        for (int i = 0; i < numSplit; ++i) {
            term = target;
            RealType norm1 = term.normInf();
            for (int n = 1; n <= numTaylorTerm; ++n) {
                term *= reciprocal(RealType(numSplit * n));
                buffer = mat * term - traceMu * term;
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
    constexpr typename MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>::TrivialType
    MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>::calcTheta(int numTaylorTerm) {
        assert(1 <= numTaylorTerm && numTaylorTerm <= MaxNumTaylorTerm && "[Error]: Invalid param");
        const int bufferIndex = (numTaylorTerm - 1) / 5; 
        return IsFloat ? theta_single[bufferIndex] : theta_double[bufferIndex];
    }

    template<class MatrixType, class VectorType>
    std::pair<int, int> MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>::calcParam(RealType traceMu) const {
        constexpr static TrivialType NormLimit = ((2 * MaxNormOrder * (MaxNormOrder + 3)) * calcTheta(MaxNumTaylorTerm)) / MaxNumTaylorTerm;
        const TrivialType norm1 = (mexp.getMatrix() - UnitMatrix<ScalarType>(getLength()) * traceMu).norm1_power(MaxNormIteration);
        const bool isSmallNorm = norm1 <= NormLimit;
        int cost = std::numeric_limits<int>::max();
        int numMinCostTerm = 0;
        if (isSmallNorm) {
            int numSplit = 1;
            for (int numTerm = 1; numTerm <= MaxNumTaylorTerm; ++numTerm) {
                const int split = int(norm1 / calcTheta(numTerm)) + 1;
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
        for (int order = 2; order <= MaxNormOrder + 1; ++order)
            powerNorms[order - 2] = calcPowerNorm(traceMu, order);

        for (int order = 2; order <= MaxNormOrder; ++order) {
            const TrivialType powerNorm = std::max(powerNorms[order - 2], powerNorms[order - 1]);
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
    inline typename MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>::RealType
    MatrixVectorProduct<MatrixExp<MatrixType>, VectorType>::calcPowerNorm(RealType traceMu, int order) const {
        const RealType norm1 = pow(mexp.getMatrix() - traceMu * UnitMatrix<ScalarType>(getLength()), order).norm1_power(MaxNormIteration);
        return pow(norm1, reciprocal(RealType(order)));
    }
}
