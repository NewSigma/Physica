/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <QtCore/qlogging.h>
#include "DNN.h"
#include "Layer.h"
#include "Node.h"

using namespace Physica::Core;

namespace Physica::AI {
    /*
     * size: length of arr nodeCount, size of layers.
     * nodeCounts: number of nodes of each layers.
     */
    DNN::DNN(int inputSize, int size, int* nodeCounts)
            : layers(std::vector<Layer*>(size)), inputs(Vector::randomVector(inputSize)), expect(getZero()), learnRate(getOne()), inputSize(inputSize) {
        for(int i = 0; i < size; ++i)
            layers[i] = new Layer(i, nodeCounts[i], this);
    }

    DNN::~DNN() {
        for(auto& layer : layers)
            delete layer;
    }

    void DNN::loadData(const Vector& loadInputs, const Numerical& loadExpect) {
        if(loadInputs.getLength() != inputSize) {
            qWarning("Insufficient data!");
            return;
        }
        inputs = loadInputs;
        expect = loadExpect;
    }

    void DNN::train() const {
        int i = 0;
        for(; i < layers.size() - 1; ++i)
            layers[i]->update();
        Numerical loss(expect);
        for(int j = 0; j < layers[i]->getSize(); ++j)
            loss -= (*layers[i])[j].calc();
        loss /= Numerical(static_cast<SignedNumericalUnit>(layers[i]->getSize()));
        for(int k = 0; k < layers[i]->getSize(); ++k)
            (*layers[i])[k].acceptedLoss = loss;
        for(; i > 1; --i)
            layers[i]->handleLoss();
    }

    Numerical DNN::predict() const {
        int i = 0;
        for(; i < layers.size() - 1; ++i)
            layers[i]->update();
        Numerical loss(expect);
        for(int j = 0; j < layers[i]->getSize(); ++j)
            loss -= (*layers[i])[j].calc();
        return loss;
    }
}
