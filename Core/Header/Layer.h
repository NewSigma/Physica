#ifndef PHYSICA_LAYER_H
#define PHYSICA_LAYER_H

#include <vector>
#include "Node.h"

class Layer {
public:
    Layer(int length, int size);
    ~Layer();

    Node* getNode(int index);
    int getLength() { return length; };
    void update(Vector* variables);
private:
    Node** nodes;
    int length;
};

#endif