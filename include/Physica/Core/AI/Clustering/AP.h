/*
 * Copyright 2022 Weibo He.
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

#include <unordered_set>
#include <forward_list>
#include <algorithm>
#include "Physica/Core/Exception/BadConvergenceException.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] Frey B. J., Dueck D. Clustering by Passing Messages Between Data Points[J]. Science. 2007 Feb 16;315(5814):972-6
     */
    template<class ScalarType>
    class AP {
        using MatrixType = Core::DenseMatrix<ScalarType, Core::MatrixOption::Row | Core::MatrixOption::Vector>;
        using SimilarMatrix = Core::DenseSymmMatrix<ScalarType>;
    private:
        MatrixType responsibility;
        MatrixType availabilities;
        ScalarType mixing;
        unsigned int numConverge;
        unsigned int maxIteration;
        std::unordered_set<size_t> exemplars{};
    public:
        AP(size_t size, ScalarType mixing, unsigned numConverge_, unsigned int maxIteration_);
        AP(const SimilarMatrix& similarity, ScalarType mixing, unsigned numConverge_, unsigned int maxIteration_);
        ~AP() = default;
        /* Operations */
        void compute(const SimilarMatrix& similarity);
        /* Getters */
        [[nodiscard]] const std::unordered_set<size_t>& getExemplars() const noexcept { return exemplars; }
        [[nodiscard]] std::forward_list<size_t> getCluster(size_t center) const;
    };

    template<class ScalarType>
    AP<ScalarType>::AP(size_t size, ScalarType mixing_, unsigned int numConverge_, unsigned int maxIteration_)
            : mixing(mixing_)
            , numConverge(numConverge_)
            , maxIteration(maxIteration_) {
        assert((mixing.isZero() || mixing.isPositive()) && mixing < ScalarType(1));
        responsibility.resize(size, size);
        availabilities.resize(size, size);
    }

    template<class ScalarType>
    AP<ScalarType>::AP(const SimilarMatrix& similarity, ScalarType mixing_, unsigned int numConverge_, unsigned int maxIteration_)
            : AP(similarity.getOrder(), mixing_, numConverge_, maxIteration_) {
        compute(similarity);
    }

    template<class ScalarType>
    void AP<ScalarType>::compute(const SimilarMatrix& similarity) {
        const size_t order = similarity.getOrder();
        const ScalarType mixing2 = ScalarType(1) - mixing;

        Core::Vector<ScalarType> buffer(order);
        exemplars.clear();
        availabilities = ScalarType(0);

        unsigned int iteration = 0;
        unsigned int converge = 0;
        while (true) {
            //Update responsibility
            for (size_t r = 0; r < order; ++r) {
                for (size_t c = 0; c < order; ++c) {
                    ScalarType temp = -std::numeric_limits<ScalarType>::max();
                    for (size_t i = 0; i < order; ++i) {
                        if (i == c)
                            continue;
                        temp = std::max(temp, availabilities(r, i) + similarity(r, i));
                    }
                    const ScalarType update = similarity(r, c) - temp;
                    responsibility(r, c) = responsibility(r, c) * mixing + update * mixing2;
                }
            }
            //Update availabilities
            for (size_t r = 0; r < order; ++r) {
                for (size_t c = 0; c < order; ++c) {
                    const bool isDiag = r == c;
                    ScalarType temp = 0;
                    for (size_t i = 0; i < order; ++i) {
                        if (i == r || i == c)
                            continue;
                        temp += std::max(ScalarType(0), responsibility(i, c));
                    }
                    const ScalarType update = isDiag ? temp : std::min(ScalarType(0), responsibility(c, c) + temp);
                    availabilities(r, c) = availabilities(r, c) * mixing + update * mixing2;
                }
            }
            //Find exemplars
            for (size_t i = 0; i < order; ++i) {
                auto find_pos = exemplars.find(i);
                const bool hasExemplar = find_pos != exemplars.end();

                buffer = responsibility.row(i).asVector() + availabilities.row(i).asVector();
                const auto max_pos = std::max_element(buffer.cbegin(), buffer.cend());
                const bool isExemplar = size_t(max_pos - buffer.cbegin()) == i;

                if (hasExemplar != isExemplar) {
                    if (isExemplar)
                        exemplars.insert(i);
                    else
                        exemplars.erase(find_pos);
                    converge = 0;
                }
            }

            ++converge;
            const bool isConverged = converge > numConverge;
            if (isConverged)
                break;
            ++iteration;
            if (iteration >= maxIteration)
                throw Core::BadConvergenceException("Exceed max iteration of AP");
        };
    }

    template<class ScalarType>
    std::forward_list<size_t> AP<ScalarType>::getCluster(size_t center) const {
        assert(exemplars.find(center) != exemplars.end());
        const size_t row = responsibility.getRow();
        std::forward_list<size_t> result{};
        for (size_t r = 0; r < row; ++r) {
            if (r == center)
                continue;
            const ScalarType max = (responsibility.row(r).asVector() + availabilities.row(r)).max();
            const bool belong = (responsibility(r, center) + availabilities(r, center)) == max;
            if (belong)
                result.push_front(r);
        }
        result.push_front(center);
        return result;
    }
}
