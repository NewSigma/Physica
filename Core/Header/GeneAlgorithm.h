#ifndef _Physica_C_GeneAlgorithm_H
#define _Physica_C_GeneAlgorithm_H

#include <ctime>
#include "RealNumber.h"

class GeneAlgorithm {
public:
    enum ChooseMode {
        LinearChoose,
        RandomChoose
    };
	GeneAlgorithm(int population, RealNumber* crossoverRate, RealNumber* mutationRate);
    ~GeneAlgorithm();

	void initFunction(RealNumber* func(RealNumber*), RealNumber* lower, RealNumber* upper, ChooseMode mode);
	RealNumber** getExtremalPoint();
	void setMaxGenerations(int maxGenerations);
	void setMaxTime(int maxTime);
	void print();

private:
	//Algorithm basic args.
	//The size of initial points we should choose.
	int population;
    RealNumber* crossoverRate;
    RealNumber* mutationRate;
	//Function args.
    RealNumber* (*fitnessFunction)(RealNumber*);
    RealNumber* lowerBound;
    RealNumber* upperBound;
    RealNumber* regionLength;
	//Save
    RealNumber** points;
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