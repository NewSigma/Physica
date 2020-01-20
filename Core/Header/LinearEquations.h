#ifndef _Physica_C_LinearEquations_H
#define _Physica_C_LinearEquations_H

class LinearEquations
{
private:
	int rank;
	double* coefficient;
	double* constant;

public:
	LinearEquations(int rank, double* coefficient, double* constant);
	double* SolveCrammer();
};

#endif