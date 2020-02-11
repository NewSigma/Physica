#ifndef _Physica_C_Matrix_CUH
#define _Physica_C_Matrix_CUH

#include "device_launch_parameters.h"

class CudaMatrix
{
public:
	double* matrix;
	int row, column;

	CudaMatrix(int r, int c);
};

CudaMatrix* multiply(const CudaMatrix* m1, const CudaMatrix* m2, bool inner);
__global__ void threadMultiply(double* result, const double* m1, const double* m2, int resultColumn, long elements, int count);

CudaMatrix* add(const CudaMatrix* m1, const CudaMatrix* m2, bool inner);
__global__ void threadAdd(double* result, const double* m1, const double* m2, long elements);

#endif