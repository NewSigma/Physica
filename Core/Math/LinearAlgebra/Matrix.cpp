#include "Core/Header/Matrix.h"
#include "Core/Header/Vector.h"

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

Vector* Matrix::operator[](int i) const {
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
        auto arr = new AbstractNum*[row];
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

Matrix* Matrix::operator*(AbstractNum& n) {
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

void Matrix::operator*=(AbstractNum& n) {
    *this << *(*this * n);
}

AbstractNum* Matrix::get(int c, int r) {
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
        auto arr = new AbstractNum*[row];
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
            auto arr = new AbstractNum*[col];
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
////////////////////////////////////////Elementary Functions////////////////////////////////////////////
Matrix* reciprocal(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = reciprocal(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* sqrt(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = sqrt(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* factorial(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = factorial(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* ln(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = ln(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* log(const Matrix& n, const AbstractNum& a) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = log(*n[i], a);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* exp(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = exp(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* pow(const Matrix& n, const AbstractNum& a) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = pow(*n[i], a);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* cos(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = cos(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* sin(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = sin(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* tan(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = tan(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* sec(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = sec(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* csc(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = csc(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* cot(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = cot(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* arccos(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arccos(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* arcsin(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arcsin(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* arctan(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arctan(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* arcsec(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arcsec(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* arccsc(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arccsc(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* arccot(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arccot(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* cosh(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = cosh(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* sinh(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = sinh(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* tanh(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = tanh(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* sech(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = sech(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* csch(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = csch(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* coth(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = coth(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* arccosh(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arccosh(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* arcsinh(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arcsinh(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* arctanh(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arctanh(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* arcsech(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arcsech(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* arccsch(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arccsch(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}

Matrix* arccoth(const Matrix& n) {
    auto arr = new Vector*[n.getLength()];
    for(int i = 0; i < n.getLength(); ++i)
        arr[i] = arccoth(*n[i]);
    return new Matrix(arr, n.getLength(), n.getType());
}