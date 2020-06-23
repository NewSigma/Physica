/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Physica/AI/DNN.h"
#include "Physica/AI/Layer.h"
#include "Physica/AI/Node.h"

using namespace Physica::Core;

namespace Physica::AI {
    Layer::Layer(int id, int nodeCount, DNN* parent) : nodes(std::vector<Node*>(nodeCount)), parentNet(parent), id(id) {
        for(int i = 0; i < nodeCount; ++i)
            nodes[i] = new Node(i, this);
    }

    Layer::~Layer() {
        for(auto& node : nodes)
            delete node;
    }

    void Layer::update() {
        if(id < parentNet->getSize() - 1) {
            Layer& nextLayer = (*parentNet)[id + 1];
            for(auto node : nodes) {
                Scalar result = node->calc();
                auto connections = node->getForwardConnections();
                for(auto connection : connections)
                    nextLayer[connection.first].vector[connection.second] = result;
            }
        }
    }

    void Layer::handleLoss() {
        if(id != 0) {
            Layer& lastLayer = (*parentNet)[id - 1];
            for(int i = 0; i < lastLayer.getSize(); ++i)
                lastLayer[i].acceptedLoss = BasicConst::getInstance().get_0();
        }
        for(int j = 0; j < getSize(); ++j)
            nodes[j]->handleLoss();
    }
}