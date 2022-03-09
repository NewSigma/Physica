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

#include <torch/torch.h>
#include "RegressionDataset.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Statistics/NumCharacter.h"

namespace Physica::AI {
    template<class Dataset>
    class KFold {
        using ScalarType = Core::Scalar<Core::Double, false>;
        using Metrics = Core::DenseMatrix<ScalarType, Core::DenseMatrixOption::Row | Core::DenseMatrixOption::Vector, 2>;
        using TrainSet = Dataset;
        using ValidSet = Dataset;
        
        Dataset set;
        unsigned int numFold;
        Metrics metrics;
    public:
        KFold(Dataset set_, unsigned int numFold);
        ~KFold() = default;
        /* Operations */
        template<class Model>
        void validation(Model model);
        /* Getters */
        [[nodiscard]] ScalarType getTrainLoss() const { return mean(metrics.row(0)); }
        [[nodiscard]] ScalarType getValidLoss() const { return mean(metrics.row(1)); }
    private:
        std::pair<TrainSet, ValidSet> cutDataset(unsigned int fold);
    };

    template<class Dataset>
    KFold<Dataset>::KFold(Dataset set_, unsigned int numFold_)
            : set(std::move(set_)), numFold(numFold_), metrics(2, numFold_) {}

    template<class Dataset>
    template<class Model>
    void KFold<Dataset>::validation(Model model) {
        for (size_t i = 0; i < numFold; ++i) {
            const auto splitted_set = cutDataset(i);
            model.init();
            model.train(splitted_set.first);
            const double train_loss = model.loss(splitted_set.first.getFeatures(), splitted_set.first.getLabels());
            metrics(0, i) = ScalarType(train_loss);
            const double valid_loss = model.loss(splitted_set.second.getFeatures(), splitted_set.second.getLabels());
            metrics(1, i) = ScalarType(valid_loss);
        }
    }

    template<class Dataset>
    std::pair<typename KFold<Dataset>::TrainSet, typename KFold<Dataset>::ValidSet>
    KFold<Dataset>::cutDataset(unsigned int fold) {
        const int64_t fold_size = set.getFeatures().sizes()[0] / numFold;
        auto valid_features = set.getFeatures().slice(0, fold * fold_size, (fold + 1) * fold_size);
        auto valid_labels = set.getLabels().slice(0, fold * fold_size, (fold + 1) * fold_size);
        auto train_features1 = set.getFeatures().slice(0, 0, fold * fold_size);
        auto train_labels1 = set.getLabels().slice(0, 0, fold * fold_size);
        auto train_features2 = set.getFeatures().slice(0, (fold + 1) * fold_size, set.size().value());
        auto train_labels2 = set.getLabels().slice(0, (fold + 1) * fold_size, set.size().value());
        auto train_features = torch::cat({std::move(train_features1), std::move(train_features2)});
        auto train_labels = torch::cat({std::move(train_labels1), std::move(train_labels2)});
        return std::make_pair(TrainSet(std::move(train_features), std::move(train_labels)),
                              ValidSet(std::move(valid_features), std::move(valid_labels)));
    }
}
