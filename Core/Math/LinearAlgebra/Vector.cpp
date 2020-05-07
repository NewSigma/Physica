#include <QtCore/qlogging.h>
#include "Vector.h"
#include "Numerical.h"
Vector::Vector() : numbers(nullptr), length(0) {}
//Init a vector with random values.
Vector::Vector(int length) : numbers(new Numerical*[length]), length(length) {
    for(int i = 0; i < length; ++i)
        numbers[i] = new Numerical(randomNumerical());
}
//Convenience function to create a 2D Vector.
Vector::Vector(const Numerical& n1, const Numerical& n2) : length(2) {
    numbers = new Numerical*[2];
    numbers[0] = new Numerical(n1);
    numbers[1] = new Numerical(n2);
}
//Convenience function to create a 3D Vector.
Vector::Vector(const Numerical& n1, const Numerical& n2, const Numerical& n3) : length(3) {
    numbers = new Numerical*[3];
    numbers[0] = new Numerical(n1);
    numbers[1] = new Numerical(n2);
    numbers[2] = new Numerical(n3);
}

Vector::Vector(Numerical** n, int l) : numbers(n), length(l) {}
//Create a Vector whose angular between x-axis and itself is arg.
Vector::Vector(const Numerical& arg) : Vector(cos(arg), sin(arg)) {}

Vector::Vector(const Vector& NumericalVector) : numbers(new Numerical*[NumericalVector.length]), length(NumericalVector.length)  {}

Vector::Vector(Vector&& NumericalVector) noexcept : numbers(NumericalVector.numbers), length(NumericalVector.length)  {
    NumericalVector.numbers = nullptr;
}

Vector::Vector(const Vector* vector) : Vector(*vector) {}

Vector::~Vector() {
    if(numbers != nullptr) {
        for(int i = 0; i < length; ++i)
            delete numbers[i];
        delete[] numbers;
    }
}

Numerical Vector::toNorm() const {
    Numerical norm = getZero();
    for(int i = 0; i < length; ++i)
        norm += *numbers[i] * *numbers[i];
    return sqrt(norm);
}

void Vector::toUnit() {
    Numerical norm = toNorm();
    for(int i = 0; i < length; ++i)
        *numbers[i] = *numbers[i] / norm;
}

std::ostream& operator<<(std::ostream& os, const Vector& n) {
    os << '(';
    for(int i = 0; i < n.length - 1; ++i) {
        os << n[i] << ", ";
    }
    os << n[n.length - 1] << ')';
    return os;
}

Numerical& Vector::operator[](int i) const {
    return *numbers[i];
}

Vector& Vector::operator=(const Vector& v) noexcept {
    if(this == &v)
        return *this;
    this->~Vector();
    length = v.length;
    numbers = new Numerical*[length];
    for(int i = 0; i < length; ++i)
        numbers[i] = new Numerical(v[i]);
    return *this;
}

Vector& Vector::operator=(Vector&& v) noexcept {
    numbers = v.numbers;
    length = v.length;
    v.numbers = nullptr;
    return *this;
}

Vector Vector::operator+(const Vector& v) const {
    const Vector* longer;
    int longer_len;
    int shorter_len;
    if(length > v.length) {
        longer = this;
        longer_len = length;
        shorter_len = v.length;
    }
    else {
        longer = &v;
        longer_len = v.length;
        shorter_len = length;
    }

    auto arr = new Numerical*[longer_len];
    for(int i = 0; i < shorter_len; ++i)
        arr[i] = new Numerical(*numbers[i] + *v.numbers[i]);
    for(int j = shorter_len; j < longer_len; ++j)
        arr[j] = new Numerical(longer->numbers[j]);
    return Vector(arr, longer_len);
}

Vector Vector::operator-(const Vector& v) const {
    const Vector* longer;
    int longer_len;
    int shorter_len;
    if(length > v.length) {
        longer = this;
        longer_len = length;
        shorter_len = v.length;
    }
    else {
        longer = &v;
        longer_len = v.length;
        shorter_len = length;
    }

    auto arr = new Numerical*[longer_len];
    for(int i = 0; i < shorter_len; ++i)
        arr[i] = new Numerical(*numbers[i] - *v.numbers[i]);
    for(int j = shorter_len; j < longer_len; ++j)
        arr[j] = new Numerical(longer->numbers[j]);
    return Vector(arr, longer_len);
}

Vector Vector::operator*(const Numerical& n) const {
    auto arr = new Numerical*[length];
    for(int i = 0; i < length; ++i)
        arr[i] = new Numerical(*numbers[i] * n);
    return Vector(arr, length);
}

Numerical Vector::operator*(const Vector& v) const {
    Numerical result = getZero();
    int shorter_len = length > v.length ? v.length : length;
    for(int i = 0; i < shorter_len; ++i)
        result += *numbers[i] * *v.numbers[i];
    return result;
}
//Here the operator/ means cross product.
Vector Vector::operator/(const Vector& v) const {
    if(length == v.length) {
        if(length == 2) {
            auto arr = new Numerical*[3];
            arr[0] = new Numerical(BasicConst::getInstance().get_0());
            arr[1] = new Numerical(BasicConst::getInstance().get_0());
            arr[2] = new Numerical(*numbers[0] * *v.numbers[1] - *numbers[1] * *v.numbers[0]);
            return Vector(arr, 3);
        }

        if(length == 3) {
            auto arr = new Numerical*[3];
            arr[0] = new Numerical(*numbers[1] * *v.numbers[2] - *numbers[2] * *v.numbers[1]);
            arr[1] = new Numerical(*numbers[2] * *v.numbers[0] - *numbers[0] * *v.numbers[2]);
            arr[2] = new Numerical(*numbers[0] * *v.numbers[1] - *numbers[1] * *v.numbers[0]);
            return Vector(arr, 3);
        }
    }
    qFatal("Can not resolve the cross product of high dimensional Vector.");
}

Vector Vector::operator-() const {
    auto arr = new Numerical*[length];
    for(int i = 0; i < length; ++i)
        arr[i] = new Numerical(-*numbers[i]);
    return Vector(arr, length);
}

bool Vector::operator==(const Vector& v) const {
    if(length != v.getLength())
        return false;
    for(int i = 0; i < length; ++i)
        if(*numbers[i] != *v.numbers[i])
            return false;
    return true;
}

Numerical Vector::toArg(int axe) const {
    return toNorm() / *numbers[axe];
}
////////////////////////////////////////Elementary Functions////////////////////////////////////////////
Vector reciprocal(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(reciprocal(n[i]));
    return Vector(arr, n.getLength());
}

Vector sqrt(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(sqrt(n[i]));
    return Vector(arr, n.getLength());
}

Vector factorial(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(factorial(n[i]));
    return Vector(arr, n.getLength());
}

Vector ln(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(ln(n[i]));
    return Vector(arr, n.getLength());
}

Vector log(const Vector& n, const Numerical& a) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(log(n[i], a));
    return Vector(arr, n.getLength());
}

Vector exp(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(exp(n[i]));
    return Vector(arr, n.getLength());
}

Vector pow(const Vector& n, const Numerical& a) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(n[i] ^ a);
    return Vector(arr, n.getLength());
}

Vector cos(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(cos(n[i]));
    return Vector(arr, n.getLength());
}

Vector sin(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(sin(n[i]));
    return Vector(arr, n.getLength());
}

Vector tan(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(tan(n[i]));
    return Vector(arr, n.getLength());
}

Vector sec(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(sec(n[i]));
    return Vector(arr, n.getLength());
}

Vector csc(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(csc(n[i]));
    return Vector(arr, n.getLength());
}

Vector cot(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(cot(n[i]));
    return Vector(arr, n.getLength());
}

Vector arccos(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arccos(n[i]));
    return Vector(arr, n.getLength());
}

Vector arcsin(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arcsin(n[i]));
    return Vector(arr, n.getLength());
}

Vector arctan(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arctan(n[i]));
    return Vector(arr, n.getLength());
}

Vector arcsec(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arcsec(n[i]));
    return Vector(arr, n.getLength());
}

Vector arccsc(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arccsc(n[i]));
    return Vector(arr, n.getLength());
}

Vector arccot(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arccot(n[i]));
    return Vector(arr, n.getLength());
}

Vector cosh(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(cosh(n[i]));
    return Vector(arr, n.getLength());
}

Vector sinh(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(sinh(n[i]));
    return Vector(arr, n.getLength());
}

Vector tanh(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(tanh(n[i]));
    return Vector(arr, n.getLength());
}

Vector sech(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(sech(n[i]));
    return Vector(arr, n.getLength());
}

Vector csch(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(csch(n[i]));
    return Vector(arr, n.getLength());
}

Vector coth(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(cosh(n[i]));
    return Vector(arr, n.getLength());
}

Vector arccosh(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arccosh(n[i]));
    return Vector(arr, n.getLength());
}

Vector arcsinh(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arcsinh(n[i]));
    return Vector(arr, n.getLength());
}

Vector arctanh(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arctanh(n[i]));
    return Vector(arr, n.getLength());
}

Vector arcsech(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arcsech(n[i]));
    return Vector(arr, n.getLength());
}

Vector arccsch(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arccsch(n[i]));
    return Vector(arr, n.getLength());
}

Vector arccoth(const Vector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arccoth(n[i]));
    return Vector(arr, n.getLength());
}