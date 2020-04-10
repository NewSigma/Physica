#include "Net.h"
#include "ElementaryFunction.h"
/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
Net::Connections::Connections(Net* p_net, int s) {
    belongs = p_net;
    size = s;
    connections = new std::unordered_map<int, connection>*[size]{nullptr};
}

Net::Connections::~Connections() {
    for(int i = 0; i < size; ++i)
        delete connections[i];
    delete[] connections;
}

bool Net::Connections::connect(int from_layer, int from_node, int to_layer, int to_node, int inter) {
    void* from, *to;
    from = belongs->getLayer(from_layer);
    to = belongs->getLayer(to_layer);
    if(from == nullptr || to == nullptr)
        return false;
    from = ((Layer*)from)->getNode(from_node);
    to = ((Layer*)to)->getNode(to_node);
    if(from == nullptr || to == nullptr)
        return false;

    if(connections[from_layer] == nullptr)
        connections[from_layer] = new std::unordered_map<int, connection>();
    auto ite = connections[from_layer]->find(to_layer);
    if(ite == connections[from_layer]->end()) {
        connection connect{};
        auto inserted = connections[from_layer]->insert(std::pair<int, connection>{to_layer, connect});
        if(!inserted.second)
            return false;
        ite = inserted.first;
    }
    ite->second.insert(std::pair<int, interface>{from_layer, interface{to_layer, inter}});
    delete ((Node*)to)->vector->numbers[inter];
    ((Node*)to)->vector->numbers[inter] = ((Node*)from)->result;
    return true;
}

void Net::Connections::disconnect(int from_layer, int from_node, int to_layer, int to_node, int inter) {
    void* from, *to;
    from = belongs->getLayer(from_layer);
    to = belongs->getLayer(to_layer);
    if(from == nullptr || to == nullptr)
        return;
    from = ((Layer*)from)->getNode(from_node);
    to = ((Layer*)to)->getNode(to_node);
    if(from == nullptr || to == nullptr)
        return;

    auto ite = connections[from_layer]->find(to_layer);
    if(ite != connections[from_layer]->end()) {
        auto range = ite->second.equal_range(from_node);
        for(auto temp_ite = range.first; temp_ite != range.second; ++ite)
            if(temp_ite->second.node == to_node) {
                ite->second.erase(temp_ite);
                break;
            }
        ((Node*)to)->vector->numbers[inter] = new RealNum(randomNumerical());
    }
}
//structure[i] is number of nodes of the i th layer. This net has structure.size() layers.
Net::Net(Vector* variables, std::vector<int> structure) : variables(variables) {
    length = structure.size();
    layers = new Layer*[length];
    connections = new Connections(this, length);
    for(int i = 0; i < length; ++i)
        layers[i] = new Layer(structure[i], length);
}

Net::~Net() {
    delete connections;
    for(int i = 0; i < length; ++i)
        delete layers[i];
    delete layers;
}

Layer* Net::getLayer(int index) {
    if(index >= length)
        return nullptr;
    return layers[index];
}