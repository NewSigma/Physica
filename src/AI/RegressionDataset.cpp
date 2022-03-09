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
#include "Physica/AI/RegressionDataset.h"

namespace Physica::AI {
    RegressionDataset::RegressionDataset(torch::Tensor features_, torch::Tensor labels_)
            : features(std::move(features_)), labels(std::move(labels_)) {}
    
    c10::optional<size_t> RegressionDataset::size() const {
        const auto size_array = labels.sizes();
        if (size_array.empty())
            return c10::nullopt;
        return size_array[0];
    }

    typename RegressionDataset::ExampleType RegressionDataset::get(size_t index) {
        return ExampleType(features[index], labels[index]);
    }
}
