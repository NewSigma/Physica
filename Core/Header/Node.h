#ifndef PHYSICA_NODE_H
#define PHYSICA_NODE_H

#include "Vector.h"
#include "RealNum.h"

class Node {
public:
    Vector* vector;
    RealNum* constant;
    RealNum* result;

    Node(int size);
    ~Node();
    void update(Vector* variables);
};

#endif