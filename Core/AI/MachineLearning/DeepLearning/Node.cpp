#include "../../../Header/Node.h"
#include "../../../Header/Const.h"

extern Const_1 const_1;

Node::Node(int size) {
    auto arr = new RealNumber*[size];
    for(int i = 0; i < size; ++i)
        arr[i] = randomRealNumber();
    vector = new Vector(arr, size);
    constant = randomRealNumber();
    result = nullptr;
}

Node::~Node() {
    delete vector;
    delete constant;
    delete result;
}

void Node::update(Vector* variables) {
    delete result;
    result = *vector * *variables;
    *result += *constant;
}