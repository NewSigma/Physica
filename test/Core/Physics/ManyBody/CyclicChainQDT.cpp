/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Physics/ManyBody/DQMCImpl/CyclicChainQDT.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using RandomSource = Random<>;
constexpr double Prec = 1E-14;
constexpr int NumSite = 4;
constexpr Array<int, 5> NumSplit{2, 4, 6, 8, 16};

namespace {
    MatrixND<T> prod(const QDTDecomp<T>& qdt) {
        return qdt.getMatrixQ() * qdt.getMatrixD() * qdt.getMatrixT();
    }

    auto direct(CyclicChainQDT<T>& chain, size_t from) {
        const size_t numSplit = chain.getNumSplit();
        QDTDecomp<T> result;
        for (size_t i = 0; i < numSplit; ++i) {
            const auto& m = chain[(from + i) % numSplit];
            if (i == 0)
                result = m;
            else
                result = result * m;
        }
        return result;
    }
}

int main() {
    for (unsigned numSplit : NumSplit) {
        CyclicChainQDT<T> chain(numSplit);
        for (size_t split = 0; split < numSplit; ++split)
            chain[split] = MatrixND<T>::random_uniform<RandomSource>(NumSite);

        for (size_t split = 0; split < numSplit; ++split) {
            const size_t from = (split + 1) % numSplit;
            const size_t to = split;
            const auto& result = chain.multiply(from, to);
            const auto answer = direct(chain, from);
            expect<RandomSource>(matrixNear(prod(result), prod(answer), Prec));
            for (size_t site = 0; site < NumSite; ++site) {
                if (RandomSource::coin()) {
                    const T factor = T(1) + T::random_uniform<RandomSource>();
                    chain.single_flip(int(site), int(split), factor, reciprocal(factor));
                }
            }
            chain.invalidate(int(split));
        }
    }
    return 0;
}
