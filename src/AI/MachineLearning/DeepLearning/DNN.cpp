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
#include <QtCore/qlogging.h>
#include "Physica/AI/DNN.h"
#include "Physica/AI/Layer.h"
#include "Physica/AI/Node.h"

using namespace Physica::Core;

namespace Physica::AI {
    /*
     * size: length of arr nodeCount, size of layers.
     * nodeCounts: number of nodes of each layers.
     */
    DNN::DNN(int inputSize, int size, int* nodeCounts)
            : layers(std::vector<Layer*>(size))
            , inputs(Vector<MultiScalar>::randomVector(inputSize))
            , expect(static_cast<SignedMPUnit>(0))
            , learnRate(static_cast<SignedMPUnit>(0))
            , inputSize(inputSize) {
        for(int i = 0; i < size; ++i)
            layers[i] = new Layer(i, nodeCounts[i], this);
    }

    DNN::~DNN() {
        for(auto& layer : layers)
            delete layer;
    }

    void DNN::loadData(const Vector<MultiScalar>& loadInputs, const MultiScalar& loadExpect) {
        if(loadInputs.getLength() != static_cast<size_t>(inputSize)) {
            qWarning("Insufficient data!");
            return;
        }
        inputs = loadInputs;
        expect = loadExpect;
    }

    void DNN::train() const {
        size_t i = 0;
        for(; i < layers.size() - 1; ++i)
            layers[i]->update();
        MultiScalar loss(expect);
        for(int j = 0; j < layers[i]->getSize(); ++j)
            loss -= (*layers[i])[j].calc();
        loss /= MultiScalar(static_cast<SignedMPUnit>(layers[i]->getSize()));
        for(int k = 0; k < layers[i]->getSize(); ++k)
            (*layers[i])[k].acceptedLoss = loss;
        for(; i > 1; --i)
            layers[i]->handleLoss();
    }

    MultiScalar DNN::predict() const {
        size_t i = 0;
        for(; i < layers.size() - 1; ++i)
            layers[i]->update();
        MultiScalar loss(expect);
        for(int j = 0; j < layers[i]->getSize(); ++j)
            loss -= (*layers[i])[j].calc();
        return loss;
    }
}
