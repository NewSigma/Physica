#include "../../Header/Vector.h"
#include "../../Header/Const.h"

extern Const_1 const_1;

Vector::Vector(RealNumber** n, int l) {
    numbers = n;
    length = l;
}

Vector::~Vector() {
    for(int i = 0; i < length; ++i)
        delete numbers[i];
    delete[] numbers;
}

Vector& Vector::operator+(Vector& v) {
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
    return *new Vector(arr, longer_len);
}

Vector& Vector::operator-(Vector& v) {
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
    return *new Vector(arr, longer_len);
}

Vector& Vector::operator*(Vector& v) {
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
        arr[i] = *numbers[i] * *v.numbers[i];
    for(int j = shorter_len; j < longer_len; ++j)
        arr[j] = new RealNumber(longer->numbers[j]);
    return *new Vector(arr, longer_len);
}
//Here the operator/ means cross product.
Vector& Vector::operator/(Vector& v) {
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
            return *new Vector(arr, 3);
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
            return *new Vector(arr, 3);
        }
    }
    return *new Vector(nullptr, 0);
}