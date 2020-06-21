/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */

#ifndef PHYSICA_DNN_H
#define PHYSICA_DNN_H

#include <vector>
#include "Core/Header/Vector.h"

using Physica::Core::Scalar;
using Physica::Core::Vector;

namespace Physica::AI {
    class Node;
    class Layer;

    class DNN {
        std::vector<Layer*> layers;
        Vector inputs;
        Scalar expect;

        Scalar learnRate;
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

        void loadData(const Vector& loadInputs, const Scalar& loadExpect);
        void setLearnRate(const Scalar& n) { learnRate = n; }
        void train() const;
        [[nodiscard]] Scalar predict() const;

        friend class Node;
    };
}

#endif