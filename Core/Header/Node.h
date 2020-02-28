#ifndef PHYSICA_NODE_H
#define PHYSICA_NODE_H

#include "Vector.h"

class Node {
public:
    Vector* vector;
    RealNumber* constant;
    RealNumber* result;

    Node(int size);
    ~Node();
    void update(Vector* variables);
};

#endif