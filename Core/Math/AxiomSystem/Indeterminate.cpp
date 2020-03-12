/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "../../Header/Indeterminate.h"

Indeterminate* Indeterminate::instance = nullptr;

Indeterminate::Indeterminate() {
    type = AbstractNum::Indeterminate;
}

Indeterminate::~Indeterminate() {
    instance = nullptr;
}

Indeterminate* Indeterminate::getInstance() {
    if(instance == nullptr)
        instance = new Indeterminate();
    return instance;
}

AbstractNum* Indeterminate::operator+(const AbstractNum& n) const {
    return getInstance();
}

AbstractNum* Indeterminate::operator-(const AbstractNum& n) const {
    return getInstance();
}

AbstractNum* Indeterminate::operator*(const AbstractNum& n) const {
    return getInstance();
}

AbstractNum* Indeterminate::operator/(const AbstractNum& n) const {
    return getInstance();
}

AbstractNum* Indeterminate::operator^(const AbstractNum& n) const {
    return getInstance();
}

AbstractNum* Indeterminate::operator-() const {
    return getInstance();
}

bool Indeterminate::operator== (const AbstractNum& n) const {
    return n.getType() == AbstractNum::Indeterminate;
}