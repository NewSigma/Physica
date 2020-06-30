/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */

#ifndef PHYSICA_DNN_H
#define PHYSICA_DNN_H

#include <vector>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector.h"

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

        Layer& operator[](int i) const { return *layers[i]; }
        [[nodiscard]] const Vector<MultiScalar>& getInputs() const { return inputs; }
        [[nodiscard]] int getInputSize() const { return inputSize; }
        [[nodiscard]] int getSize() const { return layers.size(); }

        void loadData(const Vector<MultiScalar>& loadInputs, const MultiScalar& loadExpect);
        void setLearnRate(const MultiScalar& n) { learnRate = n; }
        void train() const;
        [[nodiscard]] MultiScalar predict() const;

        friend class Node;
    };
}

#endif