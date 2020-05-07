/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */

#ifndef PHYSICA_NODE_H
#define PHYSICA_NODE_H

#include <Core/Header/Vector.h>
#include <Core/Header/Numerical.h>
#include <map>
#include <set>
#include "PhysicaAI.h"

PhysicaAI_Namespace_Begin

class Layer;

class Node {
    //Contains connection between this node to the node at next layer. <index of node, index of Numerical>
    std::map<int, int> forwardConnections;
    //Contains connection between this node to the node at previous layer. <index of Numerical, index of node>
    std::map<int, int> backwardConnections;
    Vector vector;
    Numerical bias;
    Numerical acceptedLoss;
    Numerical (*activeFunc)(const Numerical&);

    Layer* parentLayer;
    int id;
public:
    explicit Node(int id, Layer* layer);
    Node(const Node&) = delete;
    Node& operator=(const Node&) = delete;

    Numerical calc() const;
    Layer* getLayer() const { return parentLayer; }
    const std::map<int, int>& getForwardConnections() { return forwardConnections; }
    const std::map<int, int>& getBackwardConnections() { return backwardConnections; }
    void setActiveFunc(Numerical (*func)(const Numerical&)) { activeFunc = func; }
    bool connect(int toNode, int toPoint);
    void disconnect(int toNode);

    void handleLoss();

    friend class DNN;
    friend class Layer;
};

PhysicaAI_Namespace_End

#endif
