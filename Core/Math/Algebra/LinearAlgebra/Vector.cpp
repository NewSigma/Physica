/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <QtCore/qlogging.h>
#include <cstring>
#include "Vector.h"
#include "Numerical.h"

namespace Physica::Core {
    Vector::Vector() : numbers(nullptr), length(0), capacity(0) {}

    Vector::Vector(size_t length) : numbers(new Numerical[length]), length(0), capacity(length) {}
    //Convenience function to create a 2D Vector.
    Vector::Vector(const Numerical& n1, const Numerical& n2)
            : numbers(new Numerical[2]{n1, n2}), length(2), capacity(2) {}
    //Convenience function to create a 3D Vector.
    Vector::Vector(const Numerical& n1, const Numerical& n2, const Numerical& n3)
            : numbers(new Numerical[3]{n1, n2, n3}), length(3), capacity(3) {}

    Vector::Vector(Numerical*& n, size_t length) : numbers(n), length(length), capacity(length) {
        n = nullptr;
    }

    Vector::Vector(Numerical*&& n, size_t length) : numbers(n), length(length), capacity(length) {}

    Vector::Vector(const Vector& vec) : Vector(vec.length)  {
        for(size_t i = 0; i < vec.length; ++i)
            numbers[i] = vec[i];
    }

    Vector::Vector(Vector&& vec) noexcept : numbers(vec.numbers), length(vec.length), capacity(vec.capacity)  {
        vec.numbers = nullptr;
    }

    Vector::~Vector() {
        delete[] numbers;
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

    void Vector::resize(size_t size) {
        auto arr = new Numerical[size];
        length = size > length ? length : size;
        for(size_t i = 0; i < length; ++i)
            arr[i] = std::move(numbers[i]);
        delete[] numbers;
        numbers = arr;
    }

    void Vector::squeeze() {
        if(length < capacity) {
            auto arr = new Numerical[length];
            for(size_t i = 0; i < length; ++i)
                arr[i] = std::move(numbers[i]);
            delete[] numbers;
            numbers = arr;
            capacity = length;
        }
    }

    void Vector::append(Numerical n) noexcept {
        const bool doAlloc = length == capacity;
        auto arr = doAlloc ? new Numerical[length + 1] : numbers;
        if(doAlloc) {
            for(size_t i = 0; i < length; ++i)
                arr[i] = std::move(numbers[i]);
            delete[] numbers;
            ++capacity;
        }
        arr[length] = std::move(n);
        numbers = arr;
        ++length;
    }

    void Vector::append(Vector v) noexcept {
        const auto new_length = length + v.length;
        const bool doAlloc = new_length > capacity;
        auto new_arr = doAlloc ? new Numerical[new_length] : numbers;
        if(doAlloc) {
            for(size_t i = 0; i < length; ++i)
                new_arr[i] = std::move(numbers[i]);
            delete[] numbers;
        }
        for(size_t i = length; i < new_length; ++i)
            new_arr[i] = std::move(v.numbers[i - length]);
        numbers = new_arr;
        length = new_length;
    }

    Vector Vector::cut(size_t from) {
        auto new_length = length - from;
        length = from;
        auto arr = new Numerical[new_length];
        for(size_t i = 0; from < length; ++from, ++i)
            arr[i] = std::move(numbers[from]);
        return Vector(arr, new_length);
    }

    Numerical Vector::cutLast() {
        --length;
        return Numerical(std::move(numbers[length - 1]));
    }
    /*
     * Return the sub vector of current vector. From is included and to is excluded.
     */
    Vector Vector::subVector(size_t from, size_t to) {
        //Discard __restrict
        Vector temp((Numerical*)numbers + from, to - from);
        Vector result(temp);
        temp.numbers = nullptr;
        return result;
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

    std::ostream& operator<<(std::ostream& os, const Vector& n) {
        os << '(';
        for(size_t i = 0; i < n.length - 1; ++i)
            os << double(n[i]) << ", ";
        os << double(n[n.length - 1]) << ')';
        return os;
    }

    Vector& Vector::operator=(const Vector& v) noexcept {
        if(this == &v)
            return *this;
        this->~Vector();
        length = v.length;
        capacity = v.capacity;
        numbers = new Numerical[capacity];
        for(size_t i = 0; i < length; ++i)
            numbers[i] = v[i];
        return *this;
    }

    Vector& Vector::operator=(Vector&& v) noexcept {
        this->~Vector();
        numbers = v.numbers;
        length = v.length;
        capacity = v.capacity;
        v.numbers = nullptr;
        return *this;
    }
    //Here the operator/ means cross product.
    Vector Vector::operator/(const Vector& v) const {
        if(length == v.length) {
            if(length == 2) {
                auto arr = new Numerical[3];
                arr[0] = BasicConst::getInstance().get_0();
                arr[1] = BasicConst::getInstance().get_0();
                arr[2] = numbers[0] * v.numbers[1] - numbers[1] * v.numbers[0];
                return Vector(arr, 3);
            }

            if(length == 3) {
                auto arr = new Numerical[3];
                arr[0] = numbers[1] * v.numbers[2] - numbers[2] * v.numbers[1];
                arr[1] = numbers[2] * v.numbers[0] - numbers[0] * v.numbers[2];
                arr[2] = numbers[0] * v.numbers[1] - numbers[1] * v.numbers[0];
                return Vector(arr, 3);
            }
        }
        qFatal("Can not resolve the cross product of high dimensional Vector.");
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

    Vector randomVector(size_t length) {
        auto arr = new Numerical[length];
        for(size_t i = 0; i < length; ++i)
            arr[i] = randomNumerical();
        return Vector(arr, length);
    }

    Vector simplyMultiply(const Vector& v1, const Vector& v2) {
        const auto length = v1.getLength();
        Vector result(length);
        for(size_t i = 0; i < length; ++i)
            result[i] = v1[i] * v2[i];
        return result;
    }

    Vector operator-(const Vector& v) {
        const auto length = v.getLength();
        auto arr = new Numerical[length];
        for(size_t i = 0; i < length; ++i)
            arr[i] = -v[i];
        return Vector(arr, length);
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
        auto arr = new Numerical[length];
        for(size_t i = 0; i < length; ++i)
            arr[i] = v1[i] + v2[i];
        return Vector(arr, length);
    }

    Vector operator-(const Vector& v1, const Vector& v2) {
        const auto length = v1.getLength();
        auto arr = new Numerical[length];
        for(size_t i = 0; i < length; ++i)
            arr[i] = v1[i] - v2[i];
        return Vector(arr, length);
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
        auto arr = new Numerical[length];
        for(size_t i = 0; i < length; ++i)
            arr[i] = v[i] + n;
        return Vector(arr, length);
    }

    Vector operator-(const Vector& v, const Numerical& n) {
        const auto length = v.getLength();
        auto arr = new Numerical[length];
        for(size_t i = 0; i < length; ++i)
            arr[i] = v[i] - n;
        return Vector(arr, length);
    }

    Vector operator*(const Vector& v, const Numerical& n) {
        const auto length = v.getLength();
        auto arr = new Numerical[length];
        for(size_t i = 0; i < length; ++i)
            arr[i] = v[i] * n;
        return Vector(arr, length);
    }

    Vector operator/(const Vector& v, const Numerical& n) {
        const auto length = v.getLength();
        auto arr = new Numerical[length];
        for(size_t i = 0; i < length; ++i)
            arr[i] = v[i] / n;
        return Vector(arr, length);
    }
////////////////////////////////////////Elementary Functions////////////////////////////////////////////
    Vector reciprocal(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = reciprocal(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector sqrt(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = sqrt(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector factorial(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = factorial(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector ln(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = ln(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector log(const Vector& n, const Numerical& a) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = log(n[i], a);
        return Vector(arr, n.getLength());
    }

    Vector exp(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = exp(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector pow(const Vector& n, const Numerical& a) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = n[i] ^ a;
        return Vector(arr, n.getLength());
    }

    Vector cos(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = cos(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector sin(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = sin(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector tan(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = tan(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector sec(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = sec(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector csc(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = csc(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector cot(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = cot(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector arccos(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = arccos(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector arcsin(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = arcsin(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector arctan(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = arctan(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector arcsec(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = arcsec(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector arccsc(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = arccsc(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector arccot(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = arccot(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector cosh(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = cosh(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector sinh(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = sinh(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector tanh(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = tanh(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector sech(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = sech(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector csch(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = csch(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector coth(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = cosh(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector arccosh(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = arccosh(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector arcsinh(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = arcsinh(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector arctanh(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = arctanh(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector arcsech(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = arcsech(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector arccsch(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = arccsch(n[i]);
        return Vector(arr, n.getLength());
    }

    Vector arccoth(const Vector& n) {
        auto arr = new Numerical[n.getLength()];
        for(size_t i = 0; i < n.getLength(); ++i)
            arr[i] = arccoth(n[i]);
        return Vector(arr, n.getLength());
    }
}