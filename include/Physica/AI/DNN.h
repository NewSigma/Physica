/*
 * Copyright 2019 WeiBo He.
 *
 * This file is part of Physica.

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

#include <vector>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

using Physica::Core::MultiScalar;
using Physica::Core::Vector;

namespace Physica::AI {
    class Node;
    class Layer;

    class DNN {
        std::vector<Layer*> layers;
        Vector<MultiScalar> inputs;
        MultiScalar expect;

        MultiScalar learnRate;
        int inputSize;
    public:
        DNN(int inputSize, int size, int* nodeCounts);
        ~DNN();
        DNN(const DNN&) = delete;
        DNN& operator=(const DNN&) = delete;
        /* Operators */
        Layer& operator[](int i) const { return *layers[i]; }
        /* Operations */
        void loadData(const Vector<MultiScalar>& loadInputs, const MultiScalar& loadExpect);
        void train() const;
        [[nodiscard]] MultiScalar predict() const;
        /* Getters */
        [[nodiscard]] const Vector<MultiScalar>& getInputs() const { return inputs; }
        [[nodiscard]] int getInputSize() const { return inputSize; }
        [[nodiscard]] int getSize() const { return layers.size(); }
        /* Setters */
        void setLearnRate(const MultiScalar& n) { learnRate = n; }
        /* Friends */
        friend class Node;
    };
}
