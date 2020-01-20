#ifndef _Physica_C_Matrix_CUH
#define _Physica_C_Matrix_CUH

#include "device_launch_parameters.h"

class Matrix
{
public:
	double* matrix;
	int* row, * column;

	Matrix(int* row, int* column);
};

Matrix* multiply(const Matrix* m1, const Matrix* m2, bool inner);
__global__ void threadMultiply(double* result, double* m1, double* m2, int resultColumn, int elements, int count);

Matrix* add(const Matrix* m1, const Matrix* m2, bool inner);
__global__ void threadAdd(double* result, double* m1, double* m2, int elements);
////////////////////////////////helper function///////////////////////////////
void initializeMatrix(Matrix* m);

#endif