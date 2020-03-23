#include "../../Header/Vector.h"
#include "../../Header/Const.h"
#include "../../Header/Indeterminate.h"
#include "../../Header/Numerical.h"

Vector::Vector() : Vector(nullptr, 0) {}
//Convenience function to create a 2D Vector.
Vector::Vector(AbstractNum* n1, AbstractNum* n2) : length(2) {
    numbers = new AbstractNum*[2];
    numbers[0] = n1;
    numbers[1] = n2;
}
//Convenience function to create a 3D Vector.
Vector::Vector(AbstractNum* n1, AbstractNum* n2, AbstractNum* n3) : length(3) {
    numbers = new AbstractNum*[3];
    numbers[0] = n1;
    numbers[1] = n2;
    numbers[2] = n3;
}

Vector::Vector(AbstractNum** n, int l) {
    numbers = n;
    length = l;
}
//Create a vector whose angular between x-axis and itself is arg.
Vector::Vector(AbstractNum* arg) : Vector(new AbstractNum*[2]{ cos(*arg), sin(*arg) }, 2) {}

Vector::Vector(Vector& vector) {
    length = vector.length;
    numbers = new AbstractNum*[length];
    for(int i = 0; i < length; ++i)
        numbers[i] = vector.numbers[i]->concretize();
}

Vector::Vector(Vector* vector) : Vector(*vector) {}

Vector::~Vector() {
    for(int i = 0; i < length; ++i)
        delete numbers[i];
    delete[] numbers;
}

AbstractNum* Vector::toNorm() const {
    AbstractNum* norm = *numbers[0] * *numbers[0];
    for(int i = 1; i < length; ++i) {
        auto temp = *numbers[i] * *numbers[i];
        auto temp1 = *norm + *temp;
        delete temp;
        delete norm;
        norm = temp1;
    }
    return norm;
}

AbstractNum* Vector::toArg(int index) const {
    if(index < length) {
        auto norm = toNorm();
        auto result_cos = *numbers[index] / *norm;
        auto result = arccos(*result_cos);
        delete norm;
        delete result_cos;
        return result;
    }
    return Indeterminate::getInstance();
}

void Vector::toUnit() {
    AbstractNum* norm = toNorm();
    for(int i = 0; i < length; ++i) {
        auto temp = *numbers[i] / *norm;
        delete numbers[i];
        numbers[i] = temp;
    }
    delete norm;
}

std::ostream& operator<<(std::ostream& os, const Vector& n) {
    os << '(';
    for(int i = 0; i < n.length - 1; ++i) {
        os << *n[i] << ',' << ' ';
    }
    os << *n[n.length - 1] << ')';
    return os;
}

AbstractNum* Vector::operator[](int i) const {
    if(i >= length)
        return nullptr;
    return numbers[i];
}
//Move operator
void Vector::operator<<(Vector& v) {
    this->~Vector();
    numbers = v.numbers;
    length = v.length;
    v.numbers = nullptr;
    v.length = 0;
    delete &v;
}

Vector* Vector::operator+(Vector& v) {
    Vector* longer;
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

    auto arr = new AbstractNum*[longer_len];
    for(int i = 0; i < shorter_len; ++i)
        arr[i] = *numbers[i] + *v.numbers[i];
    for(int j = shorter_len; j < longer_len; ++j)
        arr[j] = longer->numbers[j]->concretize();
    return new Vector(arr, longer_len);
}

Vector* Vector::operator-(Vector& v) {
    Vector* longer;
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

    auto arr = new AbstractNum*[longer_len];
    for(int i = 0; i < shorter_len; ++i)
        arr[i] = *numbers[i] - *v.numbers[i];
    for(int j = shorter_len; j < longer_len; ++j)
        arr[j] = longer->numbers[j]->concretize();
    return new Vector(arr, longer_len);
}

Vector* Vector::operator*(const AbstractNum& n) {
    auto arr = new AbstractNum*[length];
    for(int i = 0; i < length; ++i)
        arr[i] = *numbers[i] * n;
    return new Vector(arr, length);
}

AbstractNum* Vector::operator*(Vector& v) {
    int shorter_len = length > v.length ? v.length : length;
    AbstractNum* result = new RealNum(getZero());
    for(int i = 0; i < shorter_len; ++i) {
        auto temp = *numbers[i] * *v.numbers[i];
        auto copy_old = result;
        result = *temp + *copy_old;
        delete temp;
        delete copy_old;
    }
    return result;
}
//Here the operator/ means cross product.
Vector* Vector::operator/(Vector& v) {
    if(length == v.length) {
        if(length == 2) {
            auto arr = new AbstractNum*[3];
            arr[0] = new RealNum(getZero());
            arr[1] = new RealNum(getZero());
            auto temp_1 = *numbers[0] * *v.numbers[1];
            auto temp_2 = *numbers[1] * *v.numbers[0];
            arr[2] = *temp_1 - *temp_2;
            delete temp_1;
            delete temp_2;
            return new Vector(arr, 3);
        }

        if(length == 3) {
            auto arr = new AbstractNum*[3];
            auto temp_1 = *numbers[1] * *v.numbers[2];
            auto temp_2 = *numbers[2] * *v.numbers[1];
            arr[0] = *temp_1 - *temp_2;
            delete temp_1;
            delete temp_2;

            temp_1 = *numbers[2] * *v.numbers[0];
            temp_2 = *numbers[0] * *v.numbers[2];
            arr[1] = *temp_1 - *temp_2;
            delete temp_1;
            delete temp_2;

            temp_1 = *numbers[0] * *v.numbers[1];
            temp_2 = *numbers[1] * *v.numbers[0];
            arr[2] = *temp_1 - *temp_2;
            delete temp_1;
            delete temp_2;
            return new Vector(arr, 3);
        }
    }
    return nullptr;
}

void Vector::operator+=(Vector& v) {
    *this << *(*this + v);
}

void Vector::operator-=(Vector& v) {
    *this << *(*this - v);
}

void Vector::operator/=(Vector& v) {
    auto divide = *this / v;
    if(divide == nullptr) {
        this->~Vector();
        numbers = nullptr;
        length = 0;
    }
    else
        *this << *divide;
}

void Vector::operator*=(AbstractNum& n) {
    *this << *(*this * n);
}

Vector* Vector::operator-() const {
    auto arr = new AbstractNum*[length];
    for(int i = 0; i < length; ++i)
        arr[i] = -*numbers[i];
    return new Vector(arr, length);
}

bool Vector::operator==(Vector& v) {
    if(length != v.getLength())
        return false;
    for(int i = 0; i < length; ++i)
        if(*numbers[i] != *v.numbers[i])
            return false;
    return true;
}
////////////////////////////////////////Elementary Functions////////////////////////////////////////////
Vector* reciprocal(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = reciprocal(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* sqrt(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = sqrt(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* factorial(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = factorial(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* ln(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = ln(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* log(const Vector& n, const AbstractNum& a) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = log(*n[i], a);
    return new Vector(arr, n.getLength());
}

Vector* exp(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = exp(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* pow(const Vector& n, const AbstractNum& a) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = pow(*n[i], a);
    return new Vector(arr, n.getLength());
}

Vector* cos(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = cos(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* sin(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = sin(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* tan(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = tan(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* sec(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = sec(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* csc(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = csc(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* cot(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = cot(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* arccos(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arccos(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* arcsin(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arcsin(*n[i]);
    return new Vector(arr, n.getLength());
}
;
Vector* arctan(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arctan(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* arcsec(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arcsec(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* arccsc(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arccsc(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* arccot(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arccot(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* cosh(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = cosh(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* sinh(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = sinh(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* tanh(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = tanh(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* sech(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = sech(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* csch(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = csch(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* coth(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = cosh(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* arccosh(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arccosh(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* arcsinh(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arcsinh(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* arctanh(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arctanh(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* arcsech(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arcsech(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* arccsch(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arccsch(*n[i]);
    return new Vector(arr, n.getLength());
}

Vector* arccoth(const Vector& n) {
    auto arr = new AbstractNum*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arccoth(*n[i]);
    return new Vector(arr, n.getLength());
}