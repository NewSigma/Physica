#ifndef PHYSICA_NODE_H
#define PHYSICA_NODE_H

#include "RealNumber.h"

class Node {
private:
    RealNumber* a;
    RealNumber* b;
public:
    Node();
    Node(RealNumber* a, RealNumber* b);
    ~Node();
    RealNumber* calculate(RealNumber* n);
};

#endif