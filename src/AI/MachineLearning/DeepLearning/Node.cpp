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
#include "Physica/AI/Node.h"
#include "Physica/AI/Layer.h"
#include "Physica/AI/DNN.h"

using namespace Physica::Core;

namespace Physica::AI {
    Node::Node(int id, Layer* parent)
            : vector(Vector<MultiScalar>::randomVector(parent->getNet()->getInputSize()))
            , bias(randomScalar<MultiPrecision, false>())
            , acceptedLoss(static_cast<SignedMPUnit>(0))
            , activeFunc(nullptr), parentLayer(parent), id(id) {}

    MultiScalar Node::calc() const {
        if(activeFunc == nullptr)
            return parentLayer->getNet()->getInputs() * vector + bias;
        else
            return activeFunc(parentLayer->getNet()->getInputs() * vector + bias);
    }
    /*!
     * Return true if connection is successfully built or false if failed.
     */
    bool Node::connect(int toNode, int toPoint) {
        DNN& net = *parentLayer->getNet();
        if(parentLayer->getId() < net.getSize() - 1) {
            Layer& nextLayer = net[parentLayer->getId() + 1];
            Node& target = nextLayer[toNode];
            if(static_cast<size_t>(toPoint) < target.vector.getLength()) {
                forwardConnections.insert(std::pair<int, int>(toNode, toPoint));
                target.backwardConnections.insert(std::pair<int, int>(toPoint, id));
                return true;
            }
        }
        return false;
    }
    /*!
     * After this operation,
     * we ensure no connection exists between *this and the node with the index toNode at next layer.
     */
    void Node::disconnect(int toNode) {
        auto ite = forwardConnections.find(toNode);
        if(ite != forwardConnections.end()) {
            DNN& net = *parentLayer->getNet();
            Layer& nextLayer = net[parentLayer->getId() + 1];
            Node& target = nextLayer[toNode];
            target.backwardConnections.erase(ite->second);
            forwardConnections.erase(ite);
        }
    }

    void Node::handleLoss() {
        MultiScalar averageLoss = acceptedLoss / MultiScalar(static_cast<SignedMPUnit>(vector.getLength()));
        for(size_t i = 0; i < vector.getLength(); ++i) {
            DNN& net = *parentLayer->getNet();
            if(backwardConnections.find(i) != backwardConnections.end()) {
                Layer& nextLayer = net[parentLayer->getId() - 1];
                Node& target = nextLayer[i];
                target.acceptedLoss += averageLoss;
                continue;
            }
            vector[i] += BasicConst::getInstance()._2 * net.learnRate * acceptedLoss * net.getInputs()[i];
        }
    }
}