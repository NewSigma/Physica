#ifndef PHYSICA_MATRIX_H
#define PHYSICA_MATRIX_H

class RealNumber;
class Vector;

class Matrix {
private:
    enum Type {
        ROW,
        COLUMN
    };

    Vector** vectors;
    int col, row;
    //When type = COLUMN, vectors stores column vectors.
    Type type;
public:
    Matrix();
    Matrix(Vector** vectors, int length, Type type = COLUMN);
    ~Matrix();

    Vector* operator[](int index);
    Matrix* operator+(Matrix& m);
    Matrix* operator-(Matrix& m);
    Matrix* operator*(Matrix& m);
    Matrix* operator/(Matrix& m);

    int getRow() { return row; };
    int getCol() { return col; };
    int getLength() { if(type == COLUMN) return row; else return col; };
    bool isEmpty() { return row == 0; };
    RealNumber* get(int c, int r);
    void toColMatrix();
    void toRowMatrix();
    void changeType();
    Matrix* inverse();
    void transpose();
};

#endif
