#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include "../../Header/Matrix.cuh"
#include "../../../Utils/Header/CudaUtil.cuh"

Matrix::Matrix(int r, int c)
{
    matrix = new double[r * c]{0};
    row = r;
    column = c;
}

//TODO 1.Set blocks and threads automatically.
///////////////////////////////////multiply////////////////////////////////////////////
/*
\param inner: when is true, return a Matrix whose matrix is a pointer that points to GPU memory so that we can use it
                       to do following calculate in GPU and avoid some memory copy.
*/
Matrix* multiply(const Matrix* m1, const Matrix* m2, bool inner) {
    Matrix* result;
    if (m1->column != m2->row) {
        printf("Error: Unmultipliable matrixes (%s: %d)", __FILE__, __LINE__);
        result = nullptr;
    }
    else {
        result = new Matrix(m1->row, m2->column);
        size_t size = m1->row * m2->column;
        size_t mem = size * sizeof(double);
        size_t mem_1 = m1->row * m1->column * sizeof(double);
        size_t mem_2 = mem_1;
        double* d_m, * d_m1, * d_m2;
        CHECK(cudaMalloc(&d_m, mem))
        CHECK(cudaMalloc(&d_m1, mem_1))
        CHECK(cudaMalloc(&d_m2, mem_2))
        CHECK(cudaMemcpy(d_m1, m1->matrix, mem_1, cudaMemcpyHostToDevice))
        CHECK(cudaMemcpy(d_m2, m2->matrix, mem_2, cudaMemcpyHostToDevice))
        dim3 block(m1->row, m2->column);
        for (int i = 0; i < m1->column; ++i) {
            threadMultiply << <1, block >> > (d_m, d_m1, d_m2, m2->column, size, i);
            CHECK(cudaDeviceSynchronize()) //useful or not?
        }

        if (inner) {
            result->matrix = d_m;
        }
        else {
            CHECK(cudaMemcpy(result->matrix, d_m, mem, cudaMemcpyDeviceToHost))
            CHECK(cudaFree(d_m))
        }
        CHECK(cudaFree(d_m1))
        CHECK(cudaFree(d_m2))
    }
    return result;
}

__global__
void threadMultiply(double* m, const double* m1, const double* m2, int resultColumn, long elements, int count) {
    if ((threadIdx.x + 1) * (threadIdx.y + 1) <= elements)
        m[threadIdx.y * resultColumn + threadIdx.x] += m1[threadIdx.y * resultColumn + count] * m2[count * resultColumn + threadIdx.x];
}
//////////////////////////////////////add/////////////////////////////////////////////
Matrix* add(const Matrix* m1, const Matrix* m2, bool inner) {
    Matrix* result;
    if (m1->row != m2->row || m1->column != m2->column) {
        printf("Error: Unaddable matrixes (%s: %d)", __FILE__, __LINE__);
        result = nullptr;
    }
    else {
        result = new Matrix(m1->row, m1->column);

        int size = m1->row * m1->column;
        size_t mem = size * sizeof(double);
        double* d_m, * d_m1, *d_m2;
        CHECK(cudaMalloc(&d_m, mem))
        CHECK(cudaMalloc(&d_m1, mem))
        CHECK(cudaMalloc(&d_m2, mem))
        CHECK(cudaMemcpy(d_m1, m1->matrix, mem, cudaMemcpyHostToDevice))
        CHECK(cudaMemcpy(d_m2, m2->matrix, mem, cudaMemcpyHostToDevice))
        dim3 block(m1->row, m1->column);
        threadAdd << <1, block>> > (d_m, d_m1, d_m2, size);
        if (inner) {
            result->matrix = d_m;
        }
        else {
            CHECK(cudaMemcpy(result->matrix, d_m, mem, cudaMemcpyDeviceToHost))
            CHECK(cudaFree(d_m))
        }
        CHECK(cudaFree(d_m1))
        CHECK(cudaFree(d_m2))
    }
    return result;
}

__global__
void threadAdd(double* result, const double* m1, const double* m2, long elements) {
    long index = threadIdx.y * blockDim.y + threadIdx.x;
    if (index < elements)
        result[index] = m1[index] + m2[index];
}