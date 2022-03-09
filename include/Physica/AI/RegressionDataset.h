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

#include <optional>
#include <torch/torch.h>

namespace Physica::AI {
    class RegressionDataset : public torch::data::datasets::Dataset<RegressionDataset> {
        torch::Tensor features;
        torch::Tensor labels;
    public:
        using Base = torch::data::datasets::Dataset<RegressionDataset>;
        using Base::ExampleType;
    public:
        RegressionDataset(torch::Tensor feature_, torch::Tensor label_);
        /* Getters */
        [[nodiscard]] c10::optional<size_t> size() const override;
        [[nodiscard]] ExampleType get(size_t index) override;
        [[nodiscard]] const torch::Tensor& getFeatures() const noexcept { return features; }
        [[nodiscard]] const torch::Tensor& getLabels() const noexcept { return labels; }
    };
}
