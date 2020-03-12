#ifndef _Physica_C_GeneAlgorithm_H
#define _Physica_C_GeneAlgorithm_H

#include <ctime>
#include "Numerical.h"

class GeneAlgorithm {
public:
    enum ChooseMode {
        LinearChoose,
        RandomChoose
    };
	GeneAlgorithm(int population, Numerical* crossoverRate, Numerical* mutationRate);
    ~GeneAlgorithm();

	void initFunction(Numerical* func(Numerical*), Numerical* lower, Numerical* upper, ChooseMode mode);
	Numerical** getExtremalPoint();
	void setMaxGenerations(int maxGenerations);
	void setMaxTime(int maxTime);
	void print();

private:
	//Algorithm basic args.
	//The size of initial points we should choose.
	int population;
    Numerical* crossoverRate;
    Numerical* mutationRate;
	//Function args.
    Numerical* (*fitnessFunction)(Numerical*);
    Numerical* lowerBound;
    Numerical* upperBound;
    Numerical* regionLength;
	//Save
    Numerical** points;
	//Stopping args
	int generations;
	clock_t startTime{};
	int maxGenerations;
	int maxTime;

	void overCross();
	void mutation();
	bool shouldStop();
};

#endif