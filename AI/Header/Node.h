/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */

#ifndef PHYSICA_NODE_H
#define PHYSICA_NODE_H

#include <Core/Header/Vector.h>
#include <Core/Header/Scalar.h>
#include <map>
#include <set>

using Physica::Core::Vector;
using Physica::Core::Scalar;

namespace Physica::AI {
    class Layer;

    class Node {
        //Contains connection between this node to the node at next layer. <index of node, index of Scalar>
        std::map<int, int> forwardConnections;
        //Contains connection between this node to the node at previous layer. <index of Scalar, index of node>
        std::map<int, int> backwardConnections;
        Vector vector;
        Scalar bias;
        Scalar acceptedLoss;
        Scalar (*activeFunc)(const Scalar&);

        Layer* parentLayer;
        int id;
    public:
        explicit Node(int id, Layer* layer);
        Node(const Node&) = delete;
        Node& operator=(const Node&) = delete;

        [[nodiscard]] Scalar calc() const;
        [[nodiscard]] Layer* getLayer() const { return parentLayer; }
        const std::map<int, int>& getForwardConnections() { return forwardConnections; }
        const std::map<int, int>& getBackwardConnections() { return backwardConnections; }
        void setActiveFunc(Scalar (*func)(const Scalar&)) { activeFunc = func; }
        bool connect(int toNode, int toPoint);
        void disconnect(int toNode);

        void handleLoss();

        friend class DNN;
        friend class Layer;
    };
}

#endif
