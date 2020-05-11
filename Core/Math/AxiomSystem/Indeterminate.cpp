/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Indeterminate.h"

namespace Physica::Core {
    Indeterminate* Indeterminate::instance = nullptr;

    Indeterminate::Indeterminate() = default;

    Indeterminate::~Indeterminate() {
        instance = nullptr;
    }

    Indeterminate* Indeterminate::getInstance() {
        if(instance == nullptr)
            instance = new Indeterminate();
        return instance;
    }

    AbstractNum::NumberType Indeterminate::getType() const noexcept {
        return AbstractNum::Indeterminate;
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
}