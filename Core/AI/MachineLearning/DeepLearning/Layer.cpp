#include "../../../Header/Layer.h"

Layer::Layer(int length, int size) : length(length) {
    nodes = new Node*[length];
    for(int i = 0; i < length; ++i)
        nodes[i] = new Node(size);
}

Layer::~Layer() {
    for(int i = 0; i < length; ++i)
        delete nodes[i];
    delete[] nodes;
}

Node* Layer::getNode(int index) {
    if(index > length - 1)
        return nullptr;
    return nodes[index];
}

void Layer::update(Vector* variables) {
    for(int i = 0; i < length; ++i) {
        nodes[i]->update(variables);
    }
}