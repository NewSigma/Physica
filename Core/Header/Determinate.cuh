#ifndef _Physica_C_Determinate_CUH
#define _Physica_C_Determinate_CUH

#include "cuda_runtime.h"

class Determinate
{
public:
	double* determinate;
	double* d_determinate;
	int rank;

	Determinate(double* determinate, int r);
	void print();
	double detGPU();
	void rowAddMultiply(int i1, int i2, double times);
	void rowNumMultiply(int i1, double times);
	void rowReplace(int i1, int i2);
	void cpyDetToDevice();
	void cpyDetToHost();
};

__global__ void rowAddMutiplyGPU(double* d, int rank, int i1, int i2, double times);
__global__ void rowNumMutiplyGPU(double* d, int rank, int i1, double times);
__global__ void rowReplaceGPU(double* d, int rank, int i1, int i2);
__global__ void rowAllReduceGPU(double* d, int rank, int i1, int index);

#endif