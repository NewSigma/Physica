#include "../../Header/Determinate.cuh"
#include "../../../Utils/Header/CudaUtil.cuh"
#include <device_launch_parameters.h>
#include <cmath>
#include <iostream>

////////////////////////////////////////////////Construction Methods////////////////////////////////////////////////
Determinate::Determinate(double* d, int r)
{
    determinate = d;
    d_determinate = nullptr;
    rank = r;
}

void Determinate::print()
{
    for (int i = 0; i < rank; i++)
    {
        //In case dumplicated calculate.
        int temp = i * rank;
        for (int j = 0; j < rank; j++)
        {
            double element = determinate[j + temp];
            //Making output more beautiful.
            if (element >= 0)
                std::cout << " ";
            std::cout << determinate[j + temp] << " ";
        }
        std::cout << "\n";
    }
}

/*
Calculate the determinate.
Disadvantage: The last element may overflow in this method.
*/
double Determinate::detGPU()
{
    cpyDetToDevice();
    for (int i = 1; i < rank; i++)
        rowAllReduceGPU << <1, rank - i >> > (d_determinate, rank, i, i);
    cpyDetToHost();

    double result = determinate[rank * rank - 1];
    for (int j = 1; j < rank - 1; j++)
    {
        result /= pow(determinate[(j - 1) * rank + j - 1], rank - 1 - j);
    }
    return result;
}

void Determinate::rowAddMultiply(int i1, int i2, double times)
{
    cpyDetToDevice();
    rowAddMutiplyGPU << <1, rank >> > (d_determinate, rank, i1, i2, times);
    cpyDetToHost();
}

void Determinate::rowNumMultiply(int i1, double times)
{
    cpyDetToDevice();
    rowNumMutiplyGPU << <1, rank >> > (d_determinate, rank, i1, times);
    cpyDetToHost();
}

void Determinate::rowReplace(int i1, int i2)
{
    cpyDetToDevice();
    rowReplaceGPU << <1, rank >> > (d_determinate, rank, i1, i2);
    cpyDetToHost();
}

void Determinate::cpyDetToDevice()
{
    unsigned int size = rank * rank * sizeof(double);
    CHECK(cudaMalloc(&d_determinate, size));
    CHECK(cudaMemcpy(d_determinate, determinate, size, cudaMemcpyHostToDevice));
}

void Determinate::cpyDetToHost()
{
    unsigned int size = rank * rank * sizeof(double);
    CHECK(cudaMemcpy(determinate, d_determinate, size, cudaMemcpyDeviceToHost));
    CHECK(cudaFree(d_determinate));
}

/////////////////////////////////////////////////Kernel Functions//////////////////////////////////////////////
/*
Mutiply row '@param i2' by '@param times' and add them to '@param i1'.
*/
__global__
void rowAddMutiplyGPU(double* det, int rank, int i1, int i2, double times)
{
    det[(i1 - 1) * rank + threadIdx.x] += det[(i2 - 1) * rank + threadIdx.x] * times;
}

/*
Mutiply a row '@param i1' by '@param times'.
*/
__global__
void rowNumMutiplyGPU(double* det, int rank, int i1, double times)
{
    det[(i1 - 1) * rank + threadIdx.x] = det[(i1 - 1) * rank + threadIdx.x] * times;
}

__global__
void rowReplaceGPU(double* det, int rank, int i1, int i2)
{
    unsigned int id1 = (i1 - 1) * rank + threadIdx.x;
    unsigned int id2 = (i2 - 1) * rank + threadIdx.x;
    double temp = det[id1];
    det[id1] = det[id2];
    det[id2] = temp;
}

/*
Reduce index.th elements in all rows that are under i1.
*/
__global__
void rowAllReduceGPU(double* det, int rank, int i1, int index)
{
    unsigned int rowToReduceId = i1 + 1 + threadIdx.x;
    double temp = det[(rowToReduceId - 1) * rank + index - 1];
    rowNumMutiplyGPU << <1, rank >> > (det, rank, rowToReduceId, det[(i1 - 1) * rank + index - 1]);
    rowAddMutiplyGPU << <1, rank >> > (det, rank, rowToReduceId, i1, -temp);
}