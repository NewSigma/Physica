#include "../../Header/LinearEquations.h"
#include <iostream>
#include <cstring>
#include "../../Header/Determinate.cuh"

LinearEquations::LinearEquations(int rank, double* coefficient, double* constant)
{
	if (rank <= 1)
	{
		std::cout << "Invalid linear equations rank.";
	}
	this->rank = rank;
	this->coefficient = coefficient;
	this->constant = constant;
}

double* LinearEquations::SolveCrammer()
{
	auto result = new double[this->rank];
	auto denominatorDeterminate = new Determinate(this->coefficient, this->rank);
	double denominator = denominatorDeterminate->detGPU();
	std::cout << denominator;
	for (int i = 0; i < this->rank; i++)
	{
		auto numeratorDeterminate = new Determinate(new double[this->rank], this->rank);
		memcpy(numeratorDeterminate->determinate, this->coefficient, (__int64_t)this->rank * this->rank * sizeof(double));
		for (int j = 0; j < this->rank; j++)
		{
			numeratorDeterminate->determinate[i + j * this->rank] = this->constant[j];
		}
		result[i] = numeratorDeterminate->detGPU() / denominator;
	}
	return result;
}