#include <QtCore/qlogging.h>
#include "Vector.h"
#include "Numerical.h"
Vector::Vector() : numbers(nullptr), length(0) {}

Vector::Vector(int length) : numbers(new Numerical[length]), length(length) {}
//Convenience function to create a 2D Vector.
Vector::Vector(const Numerical& n1, const Numerical& n2) : numbers(new Numerical[2]{n1, n2}), length(2) {}
//Convenience function to create a 3D Vector.
Vector::Vector(const Numerical& n1, const Numerical& n2, const Numerical& n3) : numbers(new Numerical[3]{n1, n2, n3}), length(3) {}

Vector::Vector(Numerical*& n, int l) : numbers(n), length(l) {
    n = nullptr;
}

Vector::Vector(Numerical*&& n, int l) : numbers(n), length(l) {}
//Create a Vector whose angular between x-axis and itself is arg.
Vector::Vector(const Numerical& arg) : Vector(cos(arg), sin(arg)) {}

Vector::Vector(const Vector& vec) : Vector(vec.length)  {
    for(int i = 0; i < vec.length; ++i)
        numbers[i] = vec[i];
}

Vector::Vector(Vector&& vec) noexcept : numbers(vec.numbers), length(vec.length)  {
    vec.numbers = nullptr;
}

Vector::Vector(const Vector* vector) : Vector(*vector) {}

Vector::~Vector() {
    delete[] numbers;
}

Numerical Vector::toNorm() const {
    Numerical norm = getZero();
    for(int i = 0; i < length; ++i)
        norm += numbers[i] * numbers[i];
    return sqrt(norm);
}

void Vector::toUnit() {
    if(isZeroVector())
        return;
    Numerical norm = toNorm();
    for(int i = 0; i < length; ++i)
        numbers[i] /= norm;
}

std::ostream& operator<<(std::ostream& os, const Vector& n) {
    os << '(';
    for(int i = 0; i < n.length - 1; ++i) {
        os << n[i] << ", ";
    }
    os << n[n.length - 1] << ')';
    return os;
}

Vector& Vector::operator=(const Vector& v) noexcept {
    if(this == &v)
        return *this;
    this->~Vector();
    length = v.length;
    numbers = new Numerical[length];
    for(int i = 0; i < length; ++i)
        numbers[i] = v[i];
    return *this;
}

Vector& Vector::operator=(Vector&& v) noexcept {
    this->~Vector();
    numbers = v.numbers;
    length = v.length;
    v.numbers = nullptr;
    return *this;
}

Vector Vector::operator+(const Vector& v) const {
    auto arr = new Numerical[length];
    for(int i = 0; i < length; ++i)
        arr[i] = numbers[i] + v.numbers[i];
    return Vector(arr, length);
}

Vector Vector::operator-(const Vector& v) const {
    auto arr = new Numerical[length];
    for(int i = 0; i < length; ++i)
        arr[i] = numbers[i] - v.numbers[i];
    return Vector(arr, length);
}

Vector Vector::operator*(const Numerical& n) const {
    auto arr = new Numerical[length];
    for(int i = 0; i < length; ++i)
        arr[i] = numbers[i] * n;
    return Vector(arr, length);
}

Numerical Vector::operator*(const Vector& v) const {
    Numerical result = getZero();
    for(int i = 0; i < length; ++i)
        result += numbers[i] * v.numbers[i];
    return result;
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

Vector Vector::operator-() const {
    auto arr = new Numerical[length];
    for(int i = 0; i < length; ++i)
        arr[i] = -numbers[i];
    return Vector(arr, length);
}

bool Vector::operator==(const Vector& v) const {
    if(length != v.getLength())
        return false;
    for(int i = 0; i < length; ++i)
        if(numbers[i] != v.numbers[i])
            return false;
    return true;
}

bool Vector::isZeroVector() const {
    if(length == 0)
        return false;
    for(int i = 0; i < length; ++i) {
        if(!numbers[i].isZero())
            return false;
    }
    return true;
}

Numerical Vector::toArg(int axe) const {
    return toNorm() / numbers[axe];
}

Vector randomVector(int length) {
    auto arr = new Numerical[length];
    for(int i = 0; i < length; ++i)
        arr[i] = randomNumerical();
    return Vector(arr, length);
}
////////////////////////////////////////Elementary Functions////////////////////////////////////////////
Vector reciprocal(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = reciprocal(n[i]);
    return Vector(arr, n.getLength());
}

Vector sqrt(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = sqrt(n[i]);
    return Vector(arr, n.getLength());
}

Vector factorial(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = factorial(n[i]);
    return Vector(arr, n.getLength());
}

Vector ln(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = ln(n[i]);
    return Vector(arr, n.getLength());
}

Vector log(const Vector& n, const Numerical& a) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = log(n[i], a);
    return Vector(arr, n.getLength());
}

Vector exp(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = exp(n[i]);
    return Vector(arr, n.getLength());
}

Vector pow(const Vector& n, const Numerical& a) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = n[i] ^ a;
    return Vector(arr, n.getLength());
}

Vector cos(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = cos(n[i]);
    return Vector(arr, n.getLength());
}

Vector sin(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = sin(n[i]);
    return Vector(arr, n.getLength());
}

Vector tan(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = tan(n[i]);
    return Vector(arr, n.getLength());
}

Vector sec(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = sec(n[i]);
    return Vector(arr, n.getLength());
}

Vector csc(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = csc(n[i]);
    return Vector(arr, n.getLength());
}

Vector cot(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = cot(n[i]);
    return Vector(arr, n.getLength());
}

Vector arccos(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arccos(n[i]);
    return Vector(arr, n.getLength());
}

Vector arcsin(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arcsin(n[i]);
    return Vector(arr, n.getLength());
}

Vector arctan(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arctan(n[i]);
    return Vector(arr, n.getLength());
}

Vector arcsec(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arcsec(n[i]);
    return Vector(arr, n.getLength());
}

Vector arccsc(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arccsc(n[i]);
    return Vector(arr, n.getLength());
}

Vector arccot(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arccot(n[i]);
    return Vector(arr, n.getLength());
}

Vector cosh(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = cosh(n[i]);
    return Vector(arr, n.getLength());
}

Vector sinh(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = sinh(n[i]);
    return Vector(arr, n.getLength());
}

Vector tanh(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = tanh(n[i]);
    return Vector(arr, n.getLength());
}

Vector sech(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = sech(n[i]);
    return Vector(arr, n.getLength());
}

Vector csch(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = csch(n[i]);
    return Vector(arr, n.getLength());
}

Vector coth(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = cosh(n[i]);
    return Vector(arr, n.getLength());
}

Vector arccosh(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arccosh(n[i]);
    return Vector(arr, n.getLength());
}

Vector arcsinh(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arcsinh(n[i]);
    return Vector(arr, n.getLength());
}

Vector arctanh(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arctanh(n[i]);
    return Vector(arr, n.getLength());
}

Vector arcsech(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arcsech(n[i]);
    return Vector(arr, n.getLength());
}

Vector arccsch(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arccsch(n[i]);
    return Vector(arr, n.getLength());
}

Vector arccoth(const Vector& n) {
    auto arr = new Numerical[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arccoth(n[i]);
    return Vector(arr, n.getLength());
}