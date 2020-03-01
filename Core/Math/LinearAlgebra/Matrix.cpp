#include "../../Header/Matrix.h"
#include "../../Header/Vector.h"

Matrix::Matrix() {
    vectors = nullptr;
    row = col = 0;
    type = COLUMN;
}

Matrix::Matrix(Vector** arr, int length, Type t) {
    vectors = arr;
    type = t;
    if(type) {
        row = vectors[0]->getLength();
        col = length;
    }
    else {
        row = length;
        col = vectors[0]->getLength();
    }
}

Matrix::Matrix(Matrix& matrix) {
    type = matrix.type;
    col = matrix.col;
    row = matrix.row;
    int length = getLength();
    vectors = new Vector*[length];
    for(int i = 0; i < length; ++i)
        vectors[i] = new Vector(*matrix.vectors[i]);
}

Matrix::~Matrix() {
    int length = getLength();
    for(int i = 0; i < length; ++i)
        delete vectors[i];
    delete[] vectors;
}

Vector* Matrix::operator[](int i) {
    if(i >= getLength())
        return nullptr;
    return vectors[i];
}

void Matrix::operator<<(Matrix& m) {
    this->~Matrix();
    vectors = m.vectors;
    col = m.col;
    row = m.row;
    type = m.type;
    m.vectors = nullptr;
    m.col = m.row = 0;
    delete &m;
}

Matrix* Matrix::operator+(Matrix& m) {
    if(row != m.row || col != m.col)
        return nullptr;
    bool thisChanged = false;
    bool mChanged = false;
    if(!type) {
        changeType();
        thisChanged = true;
    }
    if(!m.type) {
        m.changeType();
        mChanged = true;
    }

    auto new_vectors = new Vector*[col];
    for(int i = 0; i < col; ++i)
        new_vectors[i] = *vectors[i] + *m.vectors[i];

    if(thisChanged)
        changeType();
    if(mChanged)
        m.changeType();
    return new Matrix(new_vectors, row);
}

Matrix* Matrix::operator-(Matrix& m) {
    if(row != m.row || col != m.col)
        return nullptr;
    bool thisChanged = false;
    bool mChanged = false;
    if(!type) {
        changeType();
        thisChanged = true;
    }
    if(!m.type) {
        m.changeType();
        mChanged = true;
    }

    auto new_vectors = new Vector*[col];
    for(int i = 0; i < col; ++i)
        new_vectors[i] = *vectors[i] - *m.vectors[i];

    if(thisChanged)
        changeType();
    if(mChanged)
        m.changeType();
    return new Matrix(new_vectors, row);
}

Matrix* Matrix::operator*(Matrix& m) {
    if(col != m.row)
        return nullptr;
    bool thisChanged = false;
    bool mChanged = false;
    if(type) {
        changeType();
        thisChanged = true;
    }
    if(!m.type) {
        m.changeType();
        mChanged = true;
    }

    auto new_vectors = new Vector*[m.col];
    for(int i = 0; i < m.col; ++i) {
        auto arr = new RealNumber*[row];
        for(int j = 0; j < row; ++j)
            arr[j] = *vectors[j] * *m.vectors[i];
        new_vectors[i] = new Vector(arr, row);
    }

    if(thisChanged)
        changeType();
    if(mChanged)
        m.changeType();
    return new Matrix(new_vectors, row);
}
//Here the operator/ means inverse.
Matrix* Matrix::operator/(Matrix& m) {
    Matrix* inv = m.inverse();
    if(inv == nullptr)
        return nullptr;
    Matrix* result = *this * *inv;
    delete inv;
    return result;
}

Matrix* Matrix::operator*(RealNumber& n) {
    int length = getLength();
    auto arr = new Vector*[length];
    for(int i = 0; i < length; ++i)
        arr[i] = *vectors[i] * n;
    return new Matrix(arr, length, type);
}

void Matrix::operator+=(Matrix& m) {
    *this << *(*this + m);
}

void Matrix::operator-=(Matrix& m) {
    *this << *(*this - m);
}

void Matrix::operator*=(Matrix& m) {
    *this << *(*this * m);
}

void Matrix::operator/=(Matrix& m) {
    auto divide = *this / m;
    if(divide == nullptr) {
        this->~Matrix();
        vectors = nullptr;
        col = row = 0;
    }
    else
        *this << *divide;
}

void Matrix::operator*=(RealNumber& n) {
    *this << *(*this * n);
}

RealNumber* Matrix::get(int c, int r) {
    if(type)
        return (*vectors[c])[r];
    else
        return (*vectors[r])[c];
}

void Matrix::toColMatrix() {
    if(type)
        return;
    auto new_vectors = new Vector*[col];
    for(int i = 0; i < col; ++i) {
        auto arr = new RealNumber*[row];
        for(int j = 0; j < row; ++j) {
            arr[j] = vectors[j]->numbers[i];
            vectors[j]->numbers[i] = nullptr;
        }
        new_vectors[i] = new Vector(arr, row);
    }
    this->~Matrix();
    vectors = new_vectors;
    type = COLUMN;
}

void Matrix::toRowMatrix() {
    if(type) {
        auto new_vectors = new Vector*[row];
        for(int i = 0; i < row; ++i) {
            auto arr = new RealNumber*[col];
            for(int j = 0; j < col; ++j) {
                arr[j] = vectors[j]->numbers[i];
                vectors[j]->numbers[i] = nullptr;
            }
            new_vectors[i] = new Vector(arr, col);
        }
        this->~Matrix();
        vectors = new_vectors;
        type = ROW;
    }
}

void Matrix::changeType() {
    if(type)
        toRowMatrix();
    else
        toColMatrix();
}

Matrix* Matrix::inverse() {

}

void Matrix::transpose() {
    changeType();
    if(type)
        type = ROW;
    else
        type = COLUMN;

    int temp = row;
    row = col;
    col = temp;
}