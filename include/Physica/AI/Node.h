/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */

#ifndef PHYSICA_NODE_H
#define PHYSICA_NODE_H

#include <map>
#include <set>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector.h>

using Physica::Core::Vector;
using Physica::Core::MultiScalar;

namespace Physica::AI {
    class Layer;

    class Node {
        //Contains connection between this node to the node at next layer. <index of node, index of MultiScalar>
        std::map<int, int> forwardConnections;
        //Contains connection between this node to the node at previous layer. <index of MultiScalar, index of node>
        std::map<int, int> backwardConnections;
        Vector vector;
        MultiScalar bias;
        MultiScalar acceptedLoss;
        MultiScalar (*activeFunc)(const MultiScalar&);

        Layer* parentLayer;
        int id;
    public:
        explicit Node(int id, Layer* layer);
        Node(const Node&) = delete;
        Node& operator=(const Node&) = delete;

        [[nodiscard]] MultiScalar calc() const;
        [[nodiscard]] Layer* getLayer() const { return parentLayer; }
        const std::map<int, int>& getForwardConnections() { return forwardConnections; }
        const std::map<int, int>& getBackwardConnections() { return backwardConnections; }
        void setActiveFunc(MultiScalar (*func)(const MultiScalar&)) { activeFunc = func; }
        bool connect(int toNode, int toPoint);
        void disconnect(int toNode);

        void handleLoss();

        friend class DNN;
        friend class Layer;
    };
}

#endif
