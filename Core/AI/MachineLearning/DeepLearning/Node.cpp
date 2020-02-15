#include "../../../Header/Node.h"
#include "../../../Header/Const.h"

extern Const_1 const_1;

Node::Node() : Node(randomRealNumber(), randomRealNumber()) {}

Node::Node(RealNumber* a, RealNumber* b) : a(a), b(b) {}

Node::~Node() {
    delete a;
    delete b;
}

RealNumber* Node::calculate(RealNumber* n) {
    auto result = *a * *n;
    *result += *b;
    return result;
}