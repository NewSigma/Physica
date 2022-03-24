/*
 * Copyright 2022 WeiBo He.
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

#include <random>
#include "KFold.h"
#include "Physica/AI/Metrics.h"

namespace Physica::AI {
    template<class ModelType>
    class RandomSearch {
        using DataSet = typename ModelType::DataSet;
        using HyperParams = typename ModelType::HyperParams;
        constexpr static double uninitialized_score = -1;

        HyperParams params;
        double score;
    public:
        RandomSearch();
        RandomSearch(HyperParams params_, double score_);
        /* Operations */
        template<class RandomGenerator>
        void search(unsigned int step, Model<ModelType>& model, KFold<DataSet>& kFold, RandomGenerator& gen);
        /* Getters */
        [[nodiscard]] const HyperParams& getParams() const noexcept { return params; }
        [[nodiscard]] double getScore() const noexcept { return score; }
    };

    template<class ModelType>
    RandomSearch<ModelType>::RandomSearch() : params(), score(uninitialized_score) {}

    template<class ModelType>
    RandomSearch<ModelType>::RandomSearch(HyperParams params_, double score_) : params(std::move(params_)), score(score_) {}

    template<class ModelType>
    template<class RandomGenerator>
    void RandomSearch<ModelType>::search(unsigned int step, Model<ModelType>& model, KFold<DataSet>& kFold, RandomGenerator& gen) {
        for (unsigned int _ = 0; _ < step; ++_) {
            model.active_params = HyperParams::randomSet(gen);
            kFold.validate(model);
            const double new_score = mixed_loss(double(kFold.getTrainLoss()), double(kFold.getValidLoss()));
            if (new_score < score || score == uninitialized_score) {
                score = new_score;
                params = model.active_params;
            }
        }
    }
}
