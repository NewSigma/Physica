#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include "../../Header/Matrix.cuh"
#include "../../../Utils/Header/CudaUtil.cuh"

Matrix::Matrix(int* row, int* column)
{
	matrix = new double[(__int64_t)(*row) * (*column)];
	(*this).row = row;
	(*this).column = column;
}

//TODO 1.Set blocks and threads automatically.
///////////////////////////////////multiply////////////////////////////////////////////
/*
\param inner: when is true, return a Matrix whose matrix is a pointer that points to GPU memory so that we can use it
                       to do following calculate in GPU and avoid some memory copy.
*/
Matrix* multiply(const Matrix* m1, const Matrix* m2, bool inner)
{
	if (m1->column != m2->row)
	{
		printf("Error: Unmultipliable matrixes (%s: %d)", __FILE__, __LINE__);
		exit(1);
	}
	auto m = new Matrix(m1->row, m2->column);
	initializeMatrix(m);
	int size = *(m1->row) * *(m2->column),
		 mem = size * sizeof(double), mem_1 = *(m1->row) * *(m1->column) * sizeof(double), mem_2 = *(m2->row) * *(m2->column) * sizeof(double);
	double* d_m, * d_m1, * d_m2;
	CHECK(cudaMalloc(&d_m, mem));
	CHECK(cudaMalloc(&d_m1, mem_1));
	CHECK(cudaMalloc(&d_m2, mem_2));
	CHECK(cudaMemcpy(d_m1, m1->matrix, mem_1, cudaMemcpyHostToDevice))
	CHECK(cudaMemcpy(d_m2, m2->matrix, mem_2, cudaMemcpyHostToDevice));
	dim3 block(*(m1->row), *(m2->column));
	for (int i = 0; i < *(m1->column); i++)
	{
		threadMultiply << <1, block >> > (d_m, d_m1, d_m2, *(m2->column), size, i);
		CHECK(cudaDeviceSynchronize()); //useful or not?
	}
	if (inner)
	{
		m->matrix = d_m;
	}
	else
	{
		CHECK(cudaMemcpy(m->matrix, d_m, mem, cudaMemcpyDeviceToHost));
		CHECK(cudaFree(d_m));
	}
	CHECK(cudaFree(d_m1));
	CHECK(cudaFree(d_m2));
	return m;
}

__global__
void threadMultiply(double* m, double* m1, double* m2, int resultColumn, int elements, int count)
{
	if ((threadIdx.x + 1) * (threadIdx.y + 1) <= elements)
	{
		m[threadIdx.y * resultColumn + threadIdx.x] += m1[threadIdx.y * resultColumn + count] * m2[count * resultColumn + threadIdx.x];
	}
}
//////////////////////////////////////add/////////////////////////////////////////////
Matrix* add(const Matrix* m1, const Matrix* m2, bool inner) {
	if (*(m1->row) != *(m2->row) || *(m1->column) != *(m2->column))
	{
		printf("Error: Unaddable matrixes (%s: %d)", __FILE__, __LINE__);
		exit(1);
	}
	int size = *(m1->row) * *(m1->column);
	int mem = size * sizeof(double);
	auto m = new Matrix(m1->row, m1->column);
	double* d_m, * d_m1, *d_m2;
	CHECK(cudaMalloc(&d_m, mem));
	CHECK(cudaMalloc(&d_m1, mem));
	CHECK(cudaMalloc(&d_m2, mem));
	CHECK(cudaMemcpy(d_m1, m1->matrix, mem, cudaMemcpyHostToDevice));
	CHECK(cudaMemcpy(d_m2, m2->matrix, mem, cudaMemcpyHostToDevice));
	dim3 block(*(m1->row), *(m1->column));
	threadAdd <<<1, block>>> (d_m, d_m1, d_m2, size);
	if (inner)
	{
		m->matrix = d_m;
	}
	else
	{
		CHECK(cudaMemcpy(m->matrix, d_m, mem, cudaMemcpyDeviceToHost));
		CHECK(cudaFree(d_m));
	}
	CHECK(cudaFree(d_m1));
	CHECK(cudaFree(d_m2));
	return m;
}

__global__
void threadAdd(double* result, double* m1, double* m2, int elements)
{
	int index = threadIdx.y * blockDim.y + threadIdx.x;
	if (index < elements)
	{
		result[index] = m1[index] + m2[index];
	}
}
//////////////////////////helper function/////////////////////////////////
//set all elements to 0
void initializeMatrix(Matrix* m)
{
	memset(m->matrix, 0, *(m->row) * *(m->column) * sizeof(double));
}