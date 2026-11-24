/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <random>
#include "Physica/Core/Math/Statistics/VarianceAnalyzer.h"
#include <ranges>
#include <utility>
#include "Test.h"

using namespace Physica;
using T = float64;
using Analyzer = VarianceAnalyzer<T>;
using RandomSource = Random<MT19937, std::mt19937::default_seed>;

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

    template<Scalar ScalarT>
    ScalarT calcBlockError(const VarianceAnalyzer<ScalarT>& analyzer) {
        const size_t numGroup = analyzer.getNumGroup();
        VectorND<ScalarT> means(numGroup);
        for (size_t g = 0; g < numGroup; ++g)
            means[g] = analyzer[g].mean();
        const ScalarT mean = means.mean();
        return means.deviation(mean) / sqrt(ScalarT(numGroup));
    }

    void blocking() {
        /* iid */ {
            constexpr size_t numSample = 100000;
            const auto samples = VectorND<T>::random_normal<RandomSource>(numSample);
            const auto analyzer = Analyzer::blocking(samples);

            expect(analyzer.has_value());
            const T error = calcBlockError(*analyzer);
            const T expected = T(1) / sqrt(T(numSample));
            expect(error > expected * T(0.5));
            expect(error < expected * T(2));
        }
        /* ar1 */ {
            constexpr size_t numSample = 200000;
            const auto noise = VectorND<T>::random_normal<RandomSource>(numSample);
            VectorND<T> samples(numSample);
            samples[0] = noise[0];
            for (size_t i = 1; i < numSample; ++i)
                samples[i] = T(0.9) * samples[i - 1] + noise[i];

            const auto analyzer = Analyzer::blocking(samples);
            expect(analyzer.has_value());
            const T naive = samples.deviation(samples.mean()) / sqrt(T(numSample));
            const T error = calcBlockError(*analyzer);
            expect((*analyzer)[0].getLength() > 1);
            expect(error > naive * T(1.5));
            expect(analyzer->calcCorrTime() > T(0.5));
        }
        /* constant */ {
            const VectorND<T> samples(1000, T(1));
            const auto analyzer = Analyzer::blocking(samples);
            expect(analyzer.has_value());
            expect(calcBlockError(*analyzer) == T(0));
        }
        /* short */ {
            constexpr size_t numSample = 100;
            constexpr T phi = T(0.99);
            const auto noise = VectorND<T>::random_normal<RandomSource>(numSample);
            VectorND<T> samples(numSample);
            samples[0] = noise[0];
            for (size_t i = 1; i < numSample; ++i)
                samples[i] = phi * samples[i - 1] + noise[i];
            expect(!Analyzer::blocking(samples).has_value());
        }
    }
}

int main() {
    range();
    corr_convention();
    blocking();
    return 0;
}
