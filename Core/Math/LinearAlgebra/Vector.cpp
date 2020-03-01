#include "../../Header/Vector.h"
#include "../../Header/Const.h"

extern Const_1 const_1;

Vector::Vector() {
    numbers = nullptr;
    length = 0;
}

Vector::Vector(RealNumber** n, int l) {
    numbers = n;
    length = l;
}

Vector::Vector(Vector& vector) {
    length = vector.length;
    numbers = new RealNumber*[length];
    for(int i = 0; i < length; ++i)
        numbers[i] = new RealNumber(vector.numbers[i]);
}

Vector::~Vector() {
    for(int i = 0; i < length; ++i)
        delete numbers[i];
    delete[] numbers;
}

RealNumber* Vector::operator[](int i) {
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

    auto arr = new RealNumber*[longer_len];
    for(int i = 0; i < shorter_len; ++i)
        arr[i] = *numbers[i] + *v.numbers[i];
    for(int j = shorter_len; j < longer_len; ++j)
        arr[j] = new RealNumber(longer->numbers[j]);
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

    auto arr = new RealNumber*[longer_len];
    for(int i = 0; i < shorter_len; ++i)
        arr[i] = *numbers[i] - *v.numbers[i];
    for(int j = shorter_len; j < longer_len; ++j)
        arr[j] = new RealNumber(longer->numbers[j]);
    return new Vector(arr, longer_len);
}

Vector* Vector::operator*(RealNumber& n) {
    auto arr = new RealNumber*[length];
    for(int i = 0; i < length; ++i)
        arr[i] = *numbers[i] * n;
    return new Vector(arr, length);
}

RealNumber* Vector::operator*(Vector& v) {
    int shorter_len = length > v.length ? v.length : length;
    auto result = const_1.getZero();
    for(int i = 0; i < shorter_len; ++i){
        auto temp = *numbers[i] * *v.numbers[i];
        *result += *temp;
        delete temp;
    }
    return result;
}
//Here the operator/ means cross product.
Vector* Vector::operator/(Vector& v) {
    if(length == v.length) {
        if(length == 2) {
            auto arr = new RealNumber*[3];
            arr[0] = const_1.getZero();
            arr[1] = const_1.getZero();
            auto temp_1 = *numbers[0] * *v.numbers[1];
            auto temp_2 = *numbers[1] * *v.numbers[0];
            arr[2] = *temp_1 - *temp_2;
            delete temp_1;
            delete temp_2;
            return new Vector(arr, 3);
        }

        if(length == 3) {
            auto arr = new RealNumber*[3];
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

void Vector::operator*=(RealNumber& n) {
    *this << *(*this * n);
}