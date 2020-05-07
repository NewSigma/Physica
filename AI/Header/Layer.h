/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_LAYER_H
#define PHYSICA_LAYER_H

#include <vector>
#include "PhysicaAI.h"

PhysicaAI_Namespace_Begin

class DNN;

class Layer {
    std::vector<Node*> nodes;
    DNN* parentNet;
    //index of current layer in the net
    int id;
public:
    explicit Layer(int id, int nodeCount, DNN* parent);
    ~Layer();
    Layer(const Layer&) = delete;
    Layer& operator=(const Layer&) = delete;

    Node& operator[](int i) { return *nodes[i]; }
    DNN* getNet() const { return parentNet; }
    int getSize() const { return nodes.size(); }
    int getId() const { return id; }

    void update();
    void handleLoss();
};

PhysicaAI_Namespace_End

#endif
