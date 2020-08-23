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