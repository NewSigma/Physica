#ifndef _Physica_C_GeneAlgorithm_H
#define _Physica_C_GeneAlgorithm_H

#include <ctime>

class Numerical;

class GeneAlgorithm {
public:
    enum ChooseMode {
        LinearChoose,
        RandomChoose
    };
    double crossoverRate = 0.6;
    double mutationRate = 0.1;

	GeneAlgorithm(Numerical* func(const Numerical&), const Numerical* lower, const Numerical* upper, int pop = 100, ChooseMode mode = RandomChoose);
    ~GeneAlgorithm();

	Numerical** getExtremalPoint();
	void setMaxGenerations(int maxGenerations);
	void setMaxTime(int maxTime);
	void print();
private:
	//Algorithm basic args.
	//The size of initial points we should choose.
	int population;
	//Function args.
    Numerical* (*fitnessFunction)(const Numerical&);
    const Numerical* lowerBound;
    const Numerical* upperBound;
    Numerical* regionLength;
	//Save
    Numerical** points;
	//Stopping args
	int generations = 0;
	clock_t startTime{};
	int maxGenerations = 100;
	int maxTime = 600000;

	void crossover();
	void mutation();
	bool shouldStop();
};

#endif