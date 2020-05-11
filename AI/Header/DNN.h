/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */

#ifndef PHYSICA_DNN_H
#define PHYSICA_DNN_H

#include <vector>
#include "Core/Header/Vector.h"

using Physica::Core::Numerical;
using Physica::Core::Vector;

namespace Physica::AI {
    class Node;
    class Layer;

    class DNN {
        std::vector<Layer*> layers;
        Vector inputs;
        Numerical expect;

        Numerical learnRate;
        int inputSize;
    public:
        DNN(int inputSize, int size, int* nodeCounts);
        ~DNN();
        DNN(const DNN&) = delete;
        DNN& operator=(const DNN&) = delete;

        Layer& operator[](int i) const { return *layers[i]; }
        [[nodiscard]] const Vector& getInputs() const { return inputs; }
        [[nodiscard]] int getInputSize() const { return inputSize; }
        [[nodiscard]] int getSize() const { return layers.size(); }

        void loadData(const Vector& loadInputs, const Numerical& loadExpect);
        void setLearnRate(const Numerical& n) { learnRate = n; }
        void train() const;
        [[nodiscard]] Numerical predict() const;

        friend class Node;
    };
}

#endif