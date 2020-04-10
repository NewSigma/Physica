#include "Node.h"
#include "ElementaryFunction.h"

Node::Node(int size) {
    auto arr = new AbstractNum*[size];
    for(int i = 0; i < size; ++i)
        arr[i] = new RealNum(randomNumerical());
    vector = new Vector(arr, size);
    constant = new RealNum(randomNumerical());
    result = nullptr;
}

Node::~Node() {
    delete vector;
    delete constant;
    delete result;
}

void Node::update(Vector* variables) {
    delete result;
    auto temp = (RealNum*)(*vector * *variables);
    result = (RealNum*)(*temp + *constant);
}