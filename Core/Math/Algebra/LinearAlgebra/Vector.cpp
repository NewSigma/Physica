/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <QtCore/qlogging.h>
#include <cstring>
#include "Core/Header/Vector.h"
#include "Core/Header/Numerical.h"

namespace Physica::Core {
    Vector::Vector() : numbers(nullptr), length(0), capacity(0) {}

    Vector::Vector(size_t length)
            : numbers(reinterpret_cast<Numerical*>(malloc(length * sizeof(Numerical)))), length(0), capacity(length) {}
    //Convenience function to create a 2D Vector.
    Vector::Vector(const Numerical& n1, const Numerical& n2)
            : numbers(reinterpret_cast<Numerical*>(malloc(sizeof(Numerical) * 2))), length(2), capacity(2) {
        new (numbers) Numerical(n1);
        new (numbers + 1) Numerical(n2);
    }
    //Convenience function to create a 3D Vector.
    Vector::Vector(const Numerical& n1, const Numerical& n2, const Numerical& n3)
            : numbers(reinterpret_cast<Numerical*>(malloc(sizeof(Numerical) * 3))), length(3), capacity(3) {
        new (numbers) Numerical(n1);
        new (numbers + 1) Numerical(n2);
        new (numbers + 2) Numerical(n3);
    }

    Vector::Vector(Numerical*& n, size_t length) : numbers(n), length(length), capacity(length) {
        n = nullptr;
    }

    Vector::Vector(Numerical*&& n, size_t length) : numbers(n), length(length), capacity(length) {}

    Vector::Vector(const Vector& vec)
            : numbers(reinterpret_cast<Numerical*>(malloc(vec.capacity * sizeof(Numerical))))
            , length(vec.length), capacity(vec.capacity)  {
        for(size_t i = 0; i < vec.length; ++i)
            new (numbers + i) Numerical(vec.numbers[i]);
    }

    Vector::Vector(Vector&& vec) noexcept : numbers(vec.numbers), length(vec.length), capacity(vec.capacity)  {
        vec.numbers = nullptr;
        vec.length = 0;
    }

    Vector::~Vector() {
        for(size_t i = 0; i < length; ++i)
            (numbers + i)->~Numerical();
        free(numbers);
    }

    std::ostream& operator<<(std::ostream& os, const Vector& v) {
        os << '(';
        for(size_t i = 0; i < v.length - 1; ++i)
            os << double(v[i]) << ", ";
        os << double(v[v.length - 1]) << ')';
        return os;
    }

    Vector& Vector::operator=(const Vector& v) noexcept {
        if(this == &v)
            return *this;
        this->~Vector();
        length = v.length;
        capacity = v.capacity;
        numbers = reinterpret_cast<Numerical*>(malloc(capacity * sizeof(Numerical)));
        for(size_t i = 0; i < length; ++i)
            new (numbers + i) Numerical(v[i]);
        return *this;
    }

    Vector& Vector::operator=(Vector&& v) noexcept {
        this->~Vector();
        numbers = v.numbers;
        length = v.length;
        capacity = v.capacity;
        v.numbers = nullptr;
        v.length = 0;
        return *this;
    }
    //Here the operator/ means cross product.
    Vector Vector::operator/(const Vector& v) const {
        if(length == v.length) {
            if(length == 2)
                return Vector(BasicConst::getInstance().get_0(), BasicConst::getInstance().get_0()
                        , numbers[0] * v.numbers[1] - numbers[1] * v.numbers[0]);

            if(length == 3)
                return Vector(numbers[1] * v.numbers[2] - numbers[2] * v.numbers[1]
                        , numbers[2] * v.numbers[0] - numbers[0] * v.numbers[2]
                        , numbers[0] * v.numbers[1] - numbers[1] * v.numbers[0]);
        }
        qFatal("Can not resolve the cross product of high dimensional Vector.");
    }

    Numerical Vector::toNorm() const {
        Numerical norm = getZero();
        for(size_t i = 0; i < length; ++i)
            norm += square(numbers[i]);
        return sqrt(norm);
    }

    void Vector::toUnit() {
        if(isZeroVector())
            return;
        Numerical norm = toNorm();
        for(size_t i = 0; i < length; ++i)
            numbers[i] /= norm;
    }
    /*
     * Return the sub vector of current vector. From is included and to is excluded.
     */
    Vector Vector::subVector(size_t from, size_t to) {
        Q_ASSERT(from < to && to <= length);
        //Discard __restrict
        Vector temp(reinterpret_cast<Numerical*>(numbers + from), to - from);
        Vector result(temp);
        temp.numbers = nullptr;
        return result;
    }

    bool Vector::isZeroVector() const {
        if(length == 0)
            return false;
        for(size_t i = 0; i < length; ++i) {
            if(!numbers[i].isZero())
                return false;
        }
        return true;
    }

    Numerical Vector::toArg(size_t axe) const {
        return toNorm() / numbers[axe];
    }

    void Vector::resize(size_t size) {
        if(length > size) {
            for(size_t i = size; size < length; ++i)
                (numbers + i)->~Numerical();
            length = size;
        }
        numbers = reinterpret_cast<Numerical*>(realloc(numbers, size * sizeof(Numerical)));
        capacity = size;
    }

    void Vector::squeeze() {
        numbers = reinterpret_cast<Numerical*>(realloc(numbers, length * sizeof(Numerical)));
        capacity = length;
    }

    void Vector::append(Numerical n) noexcept {
        if(length == capacity) {
            ++capacity;
            numbers = reinterpret_cast<Numerical*>(realloc(numbers, capacity * sizeof(Numerical)));
        }
        grow(std::move(n));
    }

    void Vector::append(Vector v) noexcept {
        const auto new_length = length + v.length;
        if(new_length > capacity) {
            capacity = new_length;
            numbers = reinterpret_cast<Numerical*>(realloc(numbers, new_length * sizeof(Numerical)));
        }
        for(size_t i = length; i < new_length; ++i)
            new (numbers + i) Numerical(std::move(v.numbers[i - length]));
        length = new_length;
    }
    //!\from is included
    Vector Vector::cut(size_t from) {
        Q_ASSERT(from < length);
        auto new_length = length - from;
        length = from;
        auto arr = reinterpret_cast<Numerical*>(malloc(new_length * sizeof(Numerical)));
        for(size_t i = 0; from < length; ++from, ++i)
            new (arr + i) Numerical(std::move(numbers[from]));
        return Vector(arr, new_length);
    }

    Numerical Vector::cutLast() {
        Q_ASSERT(length > 0);
        --length;
        return Numerical(std::move(numbers[length - 1]));
    }

    void Vector::swap(Vector& v) noexcept {
        auto temp = numbers;
        numbers = v.numbers;
        v.numbers = temp;
        auto temp_size = length;
        length = v.length;
        v.length = temp_size;
        temp_size = capacity;
        capacity = v.capacity;
        v.capacity = temp_size;
    }
    
    Vector Vector::randomVector(size_t length) {
        Vector result(length);
        const auto numbers = result.numbers;
        for(size_t i = 0; i < length; ++i)
            new (numbers + i) Numerical(randomNumerical());
        return result;
    }

    Vector Vector::simplyMultiply(const Vector& v1, const Vector& v2) {
        const auto length = v1.getLength();
        Vector result(length);
        const auto numbers = result.numbers;
        for(size_t i = 0; i < length; ++i)
            new (numbers + i) Numerical(v1[i] * v2[i]);
        return result;
    }

    Vector operator-(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(-v[i]);
        return result;
    }

    bool operator==(const Vector& v1, const Vector& v2) {
        const auto length = v1.getLength();
        if(length != v2.getLength())
            return false;
        for(size_t i = 0; i < length; ++i)
            if(v1[i] != v2[i])
                return false;
        return true;
    }

    Vector operator+(const Vector& v1, const Vector& v2) {
        const auto length = v1.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(v1[i] + v2[i]);
        return result;
    }

    Vector operator-(const Vector& v1, const Vector& v2) {
        const auto length = v1.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(v1[i] - v2[i]);
        return result;
    }

    Numerical operator*(const Vector& v1, const Vector& v2) {
        const auto length = v1.getLength();
        Numerical result = getZero();
        for(size_t i = 0; i < length; ++i)
            result += v1[i] * v2[i];
        return result;
    }

    Vector operator+(const Vector& v, const Numerical& n) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(v[i] + n);
        return result;
    }

    Vector operator-(const Vector& v, const Numerical& n) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(v[i] - n);
        return result;
    }

    Vector operator*(const Vector& v, const Numerical& n) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(v[i] * n);
        return result;
    }

    Vector operator/(const Vector& v, const Numerical& n) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(v[i] / n);
        return result;
    }
////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    Vector reciprocal(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(reciprocal(v[i]));
        return result;
    }

    Vector sqrt(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(sqrt(v[i]));
        return result;
    }

    Vector factorial(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(factorial(v[i]));
        return result;
    }

    Vector ln(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(ln(v[i]));
        return result;
    }

    Vector log(const Vector& v, const Numerical& a) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(log(v[i], a));
        return result;
    }

    Vector exp(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(exp(v[i]));
        return result;
    }

    Vector pow(const Vector& v, const Numerical& a) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(v[i] ^ a);
        return result;
    }

    Vector cos(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(cos(v[i]));
        return result;
    }

    Vector sin(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(sin(v[i]));
        return result;
    }

    Vector tan(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(tan(v[i]));
        return result;
    }

    Vector sec(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(sec(v[i]));
        return result;
    }

    Vector csc(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(csc(v[i]));
        return result;
    }

    Vector cot(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(cot(v[i]));
        return result;
    }

    Vector arccos(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(arccos(v[i]));
        return result;
    }

    Vector arcsin(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(arcsin(v[i]));
        return result;
    }

    Vector arctan(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(arctan(v[i]));
        return result;
    }

    Vector arcsec(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(arcsec(v[i]));
        return result;
    }

    Vector arccsc(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(arccsc(v[i]));
        return result;
    }

    Vector arccot(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(arccot(v[i]));
        return result;
    }

    Vector cosh(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(cosh(v[i]));
        return result;
    }

    Vector sinh(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(sinh(v[i]));
        return result;
    }

    Vector tanh(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(tanh(v[i]));
        return result;
    }

    Vector sech(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(sech(v[i]));
        return result;
    }

    Vector csch(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(csch(v[i]));
        return result;
    }

    Vector coth(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(coth(v[i]));
        return result;
    }

    Vector arccosh(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(arccosh(v[i]));
        return result;
    }

    Vector arcsinh(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(arcsinh(v[i]));
        return result;
    }

    Vector arctanh(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(arctanh(v[i]));
        return result;
    }

    Vector arcsech(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(arcsech(v[i]));
        return result;
    }

    Vector arccsch(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(arccsch(v[i]));
        return result;
    }

    Vector arccoth(const Vector& v) {
        const auto length = v.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result.grow(arccoth(v[i]));
        return result;
    }
}