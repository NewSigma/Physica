#include <QtCore/qlogging.h>
#include "NumericalVector.h"
#include "Numerical.h"
//Init a vector with random values.
NumericalVector::NumericalVector(unsigned int l) : numbers(new Numerical*[l]), length((int)l) {
    for(int i = 0; i < l; ++i)
        numbers[i] = new Numerical(randomNumerical());
}
//Convenience function to create a 2D NumericalVector.
NumericalVector::NumericalVector(const Numerical& n1, const Numerical& n2) : length(2) {
    numbers = new Numerical*[2];
    numbers[0] = new Numerical(n1);
    numbers[1] = new Numerical(n2);
}
//Convenience function to create a 3D NumericalVector.
NumericalVector::NumericalVector(const Numerical& n1, const Numerical& n2, const Numerical& n3) : length(3) {
    numbers = new Numerical*[3];
    numbers[0] = new Numerical(n1);
    numbers[1] = new Numerical(n2);
    numbers[2] = new Numerical(n3);
}

NumericalVector::NumericalVector(Numerical** n, int l) : numbers(n), length(l) {}
//Create a NumericalVector whose angular between x-axis and itself is arg.
NumericalVector::NumericalVector(const Numerical& arg) : NumericalVector(cos(arg), sin(arg)) {}

NumericalVector::NumericalVector(const NumericalVector& NumericalVector) : numbers(new Numerical*[NumericalVector.length]), length(NumericalVector.length)  {}

NumericalVector::NumericalVector(NumericalVector&& NumericalVector) noexcept : numbers(NumericalVector.numbers), length(NumericalVector.length)  {
    NumericalVector.numbers = nullptr;
}

NumericalVector::NumericalVector(const NumericalVector* vector) : NumericalVector(*vector) {}

NumericalVector::~NumericalVector() {
    if(numbers != nullptr) {
        for(int i = 0; i < length; ++i)
            delete numbers[i];
        delete[] numbers;
    }
}

Numerical NumericalVector::toNorm() const {
    Numerical norm = getZero();
    for(int i = 0; i < length; ++i)
        norm += *numbers[i] * *numbers[i];
    return sqrt(norm);
}

void NumericalVector::toUnit() {
    Numerical norm = toNorm();
    for(int i = 0; i < length; ++i)
        *numbers[i] = *numbers[i] / norm;
}

std::ostream& operator<<(std::ostream& os, const NumericalVector& n) {
    os << '(';
    for(int i = 0; i < n.length - 1; ++i) {
        os << n[i] << ", ";
    }
    os << n[n.length - 1] << ')';
    return os;
}

Numerical& NumericalVector::operator[](int i) const {
    return *numbers[i];
}

NumericalVector& NumericalVector::operator=(const NumericalVector& v) noexcept {
    if(this == &v)
        return *this;
    this->~NumericalVector();
    length = v.length;
    numbers = new Numerical*[length];
    for(int i = 0; i < length; ++i)
        numbers[i] = new Numerical(v[i]);
    return *this;
}

NumericalVector& NumericalVector::operator=(NumericalVector&& v) noexcept {
    numbers = v.numbers;
    length = v.length;
    v.numbers = nullptr;
    return *this;
}

NumericalVector NumericalVector::operator+(const NumericalVector& v) const {
    const NumericalVector* longer;
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
    return NumericalVector(arr, longer_len);
}

NumericalVector NumericalVector::operator-(const NumericalVector& v) const {
    const NumericalVector* longer;
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
    return NumericalVector(arr, longer_len);
}

NumericalVector NumericalVector::operator*(const Numerical& n) const {
    auto arr = new Numerical*[length];
    for(int i = 0; i < length; ++i)
        arr[i] = new Numerical(*numbers[i] * n);
    return NumericalVector(arr, length);
}

Numerical NumericalVector::operator*(const NumericalVector& v) const {
    Numerical result = getZero();
    int shorter_len = length > v.length ? v.length : length;
    for(int i = 0; i < shorter_len; ++i)
        result += *numbers[i] * *v.numbers[i];
    return result;
}
//Here the operator/ means cross product.
NumericalVector NumericalVector::operator/(const NumericalVector& v) const {
    if(length == v.length) {
        if(length == 2) {
            auto arr = new Numerical*[3];
            arr[0] = new Numerical(basicConst->get_0());
            arr[1] = new Numerical(basicConst->get_0());
            arr[2] = new Numerical(*numbers[0] * *v.numbers[1] - *numbers[1] * *v.numbers[0]);
            return NumericalVector(arr, 3);
        }

        if(length == 3) {
            auto arr = new Numerical*[3];
            arr[0] = new Numerical(*numbers[1] * *v.numbers[2] - *numbers[2] * *v.numbers[1]);
            arr[1] = new Numerical(*numbers[2] * *v.numbers[0] - *numbers[0] * *v.numbers[2]);
            arr[2] = new Numerical(*numbers[0] * *v.numbers[1] - *numbers[1] * *v.numbers[0]);
            return NumericalVector(arr, 3);
        }
    }
    qFatal("Can not resolve the cross product of high dimensional NumericalVector.");
}

NumericalVector NumericalVector::operator-() const {
    auto arr = new Numerical*[length];
    for(int i = 0; i < length; ++i)
        arr[i] = new Numerical(-*numbers[i]);
    return NumericalVector(arr, length);
}

bool NumericalVector::operator==(const NumericalVector& v) const {
    if(length != v.getLength())
        return false;
    for(int i = 0; i < length; ++i)
        if(*numbers[i] != *v.numbers[i])
            return false;
    return true;
}
////////////////////////////////////////Elementary Functions////////////////////////////////////////////
NumericalVector reciprocal(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(reciprocal(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector sqrt(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(sqrt(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector factorial(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(factorial(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector ln(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(ln(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector log(const NumericalVector& n, const Numerical& a) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(log(n[i], a));
    return NumericalVector(arr, n.getLength());
}

NumericalVector exp(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(exp(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector pow(const NumericalVector& n, const Numerical& a) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(pow(n[i], a));
    return NumericalVector(arr, n.getLength());
}

NumericalVector cos(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(cos(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector sin(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(sin(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector tan(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(tan(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector sec(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(sec(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector csc(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(csc(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector cot(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(cot(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector arccos(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arccos(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector arcsin(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arcsin(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector arctan(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arctan(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector arcsec(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arcsec(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector arccsc(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arccsc(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector arccot(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arccot(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector cosh(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(cosh(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector sinh(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(sinh(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector tanh(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(tanh(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector sech(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(sech(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector csch(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(csch(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector coth(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(cosh(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector arccosh(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arccosh(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector arcsinh(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arcsinh(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector arctanh(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arctanh(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector arcsech(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arcsech(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector arccsch(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arccsch(n[i]));
    return NumericalVector(arr, n.getLength());
}

NumericalVector arccoth(const NumericalVector& n) {
    auto arr = new Numerical*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = new Numerical(arccoth(n[i]));
    return NumericalVector(arr, n.getLength());
}