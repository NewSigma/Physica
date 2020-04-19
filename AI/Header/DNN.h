/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */

#ifndef PHYSICA_DNN_H
#define PHYSICA_DNN_H

#include <vector>
#include <Core/Header/NumericalVector.h>
#include "PhysicaAI.h"
#include "Layer.h"

PhysicaAI_Namespace_Begin

class DNN {
    std::vector<Layer*> layers;
    NumericalVector inputs;
    Numerical expect;

    Numerical learnRate;
    int inputSize;
public:
    DNN(int inputSize, int size, int* nodeCounts);
    ~DNN();
    DNN(const DNN&) = delete;
    DNN& operator=(const DNN&) = delete;

    Layer& operator[](int i) const { return *layers[i]; }
    const NumericalVector& getInputs() const { return inputs; }
    int getInputSize() const { return inputSize; }
    int getSize() const { return layers.size(); }

    void loadData(const NumericalVector& loadInputs, const Numerical& loadExpect);
    void setLearnRate(const Numerical& n) { learnRate = n; }
    void train() const;
    Numerical predict() const;

    friend void Node::handleLoss();
};

PhysicaAI_Namespace_End

#endif