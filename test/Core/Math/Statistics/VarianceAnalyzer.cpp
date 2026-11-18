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
#include "Physica/Core/Math/Statistics/VarianceAnalyzer.h"
#include <ranges>
#include <utility>
#include "Test.h"

using namespace Physica;
using T = float64;
using Analyzer = VarianceAnalyzer<T>;

namespace {
    void range() noexcept {
        using It = std::ranges::iterator_t<Analyzer>;
        static_assert(std::indirectly_readable<It>);
        static_assert(std::indirectly_writable<It, std::ranges::range_value_t<Analyzer>>);
        static_assert(std::incrementable<It>);
        static_assert(std::sized_sentinel_for<It, It>);
        static_assert(std::contiguous_iterator<It>);

        static_assert(std::ranges::sized_range<Analyzer>);
        static_assert(std::ranges::random_access_range<Analyzer>);
        static_assert(std::ranges::contiguous_range<Analyzer>);
        static_assert(std::ranges::common_range<Analyzer>);
        static_assert(std::ranges::viewable_range<Analyzer>);
        static_assert(std::ranges::contiguous_range<const Analyzer>);
        static_assert(std::same_as<decltype(std::ranges::data(std::declval<const Analyzer&>())), const VectorND<T>*>);

        Analyzer analyzer(2);
        analyzer[0] = VectorND<T>{0, 0, 0, 0};
        analyzer[1] = VectorND<T>{1, 1, 1, 1};
        expect(std::ranges::size(analyzer) == analyzer.getNumGroup());
        expect(std::ranges::data(analyzer) == &analyzer[0]);
        expect(std::ranges::data(std::as_const(analyzer)) == &std::as_const(analyzer[0]));

        size_t numSample = 0;
        for (const auto& group : analyzer)
            numSample += group.getLength();
        expect(numSample == analyzer.getTotalNumSample());
    }

    void corr_convention() {
        // Test correlation time convension between \class VarianceAnalyzer and \class Correlation matches
        Analyzer analyzer(2);
        analyzer[0] = VectorND<T>{0, 0, 0, 0};
        analyzer[1] = VectorND<T>{1, 1, 1, 1};
        const T numGroup = T(analyzer.getNumGroup());
        const T numSample = T(analyzer.getTotalNumSample());
        const T expectedR = (numSample - 1) / (numGroup - 1);
        expect(scalarNear(analyzer.calcRelationCoeff(), T(1), 1E-15));
        expect(scalarNear(analyzer.calcParamR(), expectedR, 1E-15));
        expect(scalarNear(analyzer.calcCorrTime(), expectedR * T(0.5), 1E-15));
    }
}

int main() {
    range();
    corr_convention();
    return 0;
}
