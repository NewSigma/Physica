/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Node.h"
#include "Layer.h"
#include "DNN.h"
#include "Core/Header/ElementaryFunction.h"

PhysicaAI_Namespace_Begin

Node::Node(int id, Layer* parent)
: vector(parent->getNet()->getInputSize()), bias(randomNumerical())
, acceptedLoss(getZero()), activeFunc(nullptr), parentLayer(parent), id(id) {
}

Numerical Node::calc() const {
    if(activeFunc == nullptr)
        return parentLayer->getNet()->getInputs() * vector + bias;
    else
        return activeFunc(parentLayer->getNet()->getInputs() * vector + bias);
}
//Return true if connection is successfully built or false if failed.
bool Node::connect(int toNode, int toPoint) {
    DNN& net = *parentLayer->getNet();
    if(parentLayer->getId() < net.getSize() - 1) {
        Layer& nextLayer = net[parentLayer->getId() + 1];
        Node& target = nextLayer[toNode];
        if(toPoint < target.vector.getLength()) {
            forwardConnections.insert(std::pair<int, int>(toNode, toPoint));
            target.backwardConnections.insert(std::pair<int, int>(toPoint, id));
            return true;
        }
    }
    return false;
}
//After this operation, we ensure no connection exists between *this and the node with the index toNode at next layer.
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
    Numerical averageLoss = acceptedLoss / Numerical(static_cast<SignedNumericalUnit>(vector.getLength()));
    for(int i = 0; i < vector.getLength(); ++i) {
        DNN& net = *parentLayer->getNet();
        if(backwardConnections.find(i) != backwardConnections.end()) {
            Layer& nextLayer = net[parentLayer->getId() - 1];
            Node& target = nextLayer[i];
            target.acceptedLoss += averageLoss;
            continue;
        }
        vector[i] += BasicConst::getInstance().get_2() * net.learnRate * acceptedLoss * net.getInputs()[i];
    }
}

PhysicaAI_Namespace_End