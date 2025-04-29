/*
 * Copyright 2023-2024 Weibo He.
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
#include <iostream>
#include <fstream>
#include "Physica/Core/ML/DecisionTree/DecisionTree.h"

using namespace Physica;
using ScalarType = float64;
using TreeType = DecisionTree<ScalarType, DecisionTreeType::Classify>;
/*
 * Reference:
 * [1] Glass Identification Data Set; https://archive.ics.uci.edu/ml/datasets/Glass+Identification
 */
int main() {
    std::ifstream fin("glass.data");
    DenseMatrix<ScalarType> dataMat(214, 11);
    for (size_t r = 0; r < dataMat.getRow(); ++r)
        for (size_t c = 0; c < dataMat.getCol(); ++c)
            fin >> dataMat(r, c);
    Array<bool> isFeatureContinuous(9, true);

    auto features = dataMat.block(0, dataMat.getRow(), 1, 9);
    auto labels = dataMat.col(10);
    TreeType::Dataset dataset{features, labels, std::move(isFeatureContinuous)};
    auto tree = TreeType::train(dataset);
    size_t count = 0;
    for (size_t i = 0; i < labels.getLength(); ++i)
        count += tree.predict(features.row(i)) == labels[i];
    std::cout << "Training accuracy: " << ScalarType(count) / ScalarType(dataMat.getRow()) << std::endl;
    return 0;
}
