/*
 * Copyright 2026 Weibo He.
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

#include <algorithm>
#include <cassert>
#include <ranges>
#include <type_traits>
#include <utility>
#include "Physica/Core/Math/Algebra/LinearAlgebra/IndexVar/Index.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/IndexVar/IndexVar.h"
#include "Physica/Core/Utils/Container/Array.h"
#include "Physica/Core/Utils/MetaProgramming.h"
#include "Physica/Core/Utils/Range.h"

namespace Physica {
    template<class X> class Ein;

    template<class LHS, class RHS>
    class Contract {
        using This = Contract<LHS, RHS>;
        using LHS_t = std::remove_cvref_t<LHS>;
        using RHS_t = std::remove_cvref_t<RHS>;

        using T = Internal::BinaryScalarOpRtnTy<typename LHS_t::ScalarType, typename RHS_t::ScalarType>::Type;
    public:
        constexpr static int NDimLHS = LHS_t::NDim;
        constexpr static int NDimRHS = RHS_t::NDim;
        constexpr static int NumIndexIn = NDimLHS + NDimRHS;
    private:
        decay_rvalue_t<LHS> lhs;
        decay_rvalue_t<RHS> rhs;
    public:
        Contract(LHS&& lhs_, RHS&& rhs_) noexcept;
        Contract(const This&) = default;
        Contract(This&&) noexcept = default;
        ~Contract() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(instanceof<Ein> auto& target) const;
        /* Getters */
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
    };

    template<class LHS, class RHS>
    Contract<LHS, RHS>::Contract(LHS&& lhs_, RHS&& rhs_) noexcept
            : lhs(std::forward<LHS>(lhs_)), rhs(std::forward<RHS>(rhs_)) {}

    template<class LHS, class RHS>
    void Contract<LHS, RHS>::assign(instanceof<Ein> auto& target) const {
        using Target = std::remove_cvref_t<decltype(target)>;
        constexpr int NumIndexOut = Target::NDim;

        auto&& lhs = getLHS();
        auto&& rhs = getRHS();
        const auto& lvars = lhs.getVars();
        const auto& rvars = rhs.getVars();
        const auto& tvars = target.getVars();

        Array<const Var*, NumIndexIn> allVars;
        std::ranges::copy_n(lvars.begin(), NDimLHS, allVars.begin());
        std::ranges::copy_n(rvars.begin(), NDimRHS, allVars.begin() + NDimLHS);

        const auto index = std::views::iota(0, NumIndexIn);
        // first[i] is the position of the first occurrence of the variable at position i
        Array<int, NumIndexIn> first;
        Array<int, NumIndexIn> count{};
        for (int i = 0; i < NumIndexIn; ++i) {
            const auto match = std::ranges::find(allVars.begin(), allVars.begin() + i, allVars[i]);
            first[i] = match == allVars.begin() + i ? i : int(match - allVars.begin());
            count[first[i]] += 1;
        }
        assert(std::ranges::all_of(count, [](int n) { return n <= 2; })
               && "[Error]: An index cannot appear more than twice in Einstein notation");

        const auto isFree = [&](int i) noexcept {
            return i == first[i] && count[i] == 1;
        };
        [[maybe_unused]] const auto numFree = std::ranges::count_if(index, isFree);
        assert(numFree == NumIndexOut && "[Error]: Number of free indices must match output");

        // Resolve the dimension of each index, checking the consistency of repeated ones
        Index<NumIndexIn> dims{};
        for (int i = 0; i < NumIndexIn; ++i) {
            const size_t dim = i < NDimLHS ? lhs.dim(i) : rhs.dim(i - NDimLHS);
            if (i == first[i])
                dims[i] = dim;
            else
                assert(dims[first[i]] == dim && "[Error]: Dimension mismatch of a contracted index");
        }

        // Map each output index to a free (count == 1) position
        Array<int, NumIndexOut> outPos{};
        for (int t = 0; t < NumIndexOut; ++t) {
            const auto found = std::ranges::find_if(index, [&](int p) { return isFree(p) && allVars[p] == tvars[t]; });
            assert(found != index.end() && "[Error]: An output index is not a free index of the contraction");
            outPos[t] = *found;
            assert(dims[outPos[t]] == target.dim(t) && "[Error]: Output dimension mismatch");
        }

        Index<NumIndexIn> contractShape;
        for (auto&& [i, f, dim] : Physica::zip(index, first, dims))
            contractShape[i] = (i == f && count[i] == 2) ? dim : 1;

        const auto outShape = target.getShape();
        const auto sizeOut = Index<NumIndexOut>::toSize(outShape);
        const auto sizeContract = Index<NumIndexIn>::toSize(contractShape);
        Index<NumIndexIn> posValue{};
        Index<NDimLHS> lhsIndex{};
        Index<NDimRHS> rhsIndex{};
        for (size_t n = 0; n < sizeOut; ++n) {
            using OutIndexType = Target::IndexType;
            const auto outIndex = OutIndexType::toIndexND(outShape, n);
            for (int t = 0; t < NumIndexOut; ++t)
                posValue[outPos[t]] = outIndex[t];

            T sum{};
            for (size_t c = 0; c < sizeContract; ++c) {
                const Index<NumIndexIn> contractIndex = Index<NumIndexIn>::toIndexND(contractShape, c);
                for (auto&& [i, f, n] : Physica::zip(index, first, count))
                    if (i == f && n == 2)
                        posValue[i] = contractIndex[i];
                for (int p = 0; p < NDimLHS; ++p)
                    lhsIndex[p] = posValue[first[p]];
                for (int q = 0; q < NDimRHS; ++q)
                    rhsIndex[q] = posValue[first[NDimLHS + q]];
                sum = fma(lhs.calc(lhsIndex), rhs.calc(rhsIndex), sum);
            }
            target[outIndex] = sum;
        }
    }

    template<class LHS, class RHS>
    auto&& Contract<LHS, RHS>::getLHS(this auto&& self) noexcept {
        return self.lhs;
    }

    template<class LHS, class RHS>
    auto&& Contract<LHS, RHS>::getRHS(this auto&& self) noexcept {
        return self.rhs;
    }
}
