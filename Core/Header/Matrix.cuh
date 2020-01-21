#ifndef _Physica_C_Matrix_CUH
#define _Physica_C_Matrix_CUH

#include "device_launch_parameters.h"

class Matrix
{
public:
	double* matrix;
	int row, column;

	Matrix(int r, int c);
};

Matrix* multiply(const Matrix* m1, const Matrix* m2, bool inner);
__global__ void threadMultiply(double* result, const double* m1, const double* m2, int resultColumn, long elements, int count);

Matrix* add(const Matrix* m1, const Matrix* m2, bool inner);
__global__ void threadAdd(double* result, const double* m1, const double* m2, long elements);

#endif