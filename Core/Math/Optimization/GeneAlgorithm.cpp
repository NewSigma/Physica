#include "GeneAlgorithm.h"
#include <iostream>
#include "ElementaryFunction.h"
#include "Numerical.h"
/*
 * Usage: GeneAlgorithm(args)->initFunction(args)->getExtremalPoint()
 * Copyright (c) 2019 NewSigma@163.com.All rights reserved.
 */
//TODO Debug
GeneAlgorithm::GeneAlgorithm(Numerical func(const Numerical&), const Numerical& lower, const Numerical& upper, int pop, ChooseMode mode) {
    population = pop;
	points = new Numerical*[population];

    fitnessFunction = func;
    lowerBound = new Numerical(lower);
    upperBound = new Numerical(upper);
    regionLength = new Numerical(upper - lower);
    //Get abs(regionLength).
    regionLength->toAbs();

    if (mode == LinearChoose) {
        Numerical stepLength = *regionLength / Numerical(population);
        for (int i = 0; i < population; i++) {
            auto byte = (unsigned long*)malloc(sizeof(long));
            byte[0] = i;
            points[i] = new Numerical(stepLength * Numerical(byte, 1, 0) + lower);
        }
    }
    else if (mode == RandomChoose) {
        for (int i = 0; i < population; i++)
            points[i] = new Numerical(randomNumerical() * *regionLength + lower);
    }
}

GeneAlgorithm::~GeneAlgorithm() {
    delete[] points;
    delete regionLength;
    delete lowerBound;
    delete upperBound;
    delete regionLength;
}
//Get the maximum point. Multiply the function by a -1 to get the minimum.
Numerical** GeneAlgorithm::getExtremalPoint() {
	if (fitnessFunction == nullptr) {
		std::cout << "Uninitialized function.\n";
		return nullptr;
	}
	startTime = clock();
	while (!shouldStop()) {
        crossover();
		mutation();
		generations += 1;
	}
	return points;
}

void GeneAlgorithm::crossover() {
	for (int i = 0; i < population; i++) {
        auto r = randomNumerical();
		if (Numerical(crossoverRate) > r) {
			long randomIndex1 = random() % population;
			long randomIndex2 = random() % population;
            auto random1 = points[randomIndex1];
            auto random2 = points[randomIndex2];

			r = randomNumerical();
            //Whether random2 > random1 or not is not important.
            auto child = (*random2 - *random1) * r + *random1;
            auto child_y = fitnessFunction(child);
            auto y_random1 = fitnessFunction(*random1);
            auto y_random2 = fitnessFunction(*random2);
			if (child_y > y_random1)
                *random1 = child;
			else if (child_y > y_random2)
                *random2 = child;
		}
	}
}

void GeneAlgorithm::mutation() {
	if (Numerical(mutationRate) > randomNumerical()) {
        long randomIndex = random() % population;
        *points[randomIndex] = randomNumerical() * *regionLength + *lowerBound;
    }
}

bool GeneAlgorithm::shouldStop() {
	if (generations > maxGenerations) {
		std::cout << "[GeneAlgorithm] Max generation reached: " << maxGenerations << std::endl;
		return true;
	}
	else if (clock() - startTime > maxTime) {
		std::cout << "[GeneAlgorithm] Max time reached: " << maxTime << std::endl;
		return true;
	}
	return false;
}

void GeneAlgorithm::setMaxGenerations(int max) {
	this->maxGenerations = max;
}

void GeneAlgorithm::setMaxTime(int time) {
	this->maxTime = time;
}

void GeneAlgorithm::print() {
	if(fitnessFunction != nullptr) {
        for (int i = 0; i < population; ++i)
            std::cout << *points[i] << " " << fitnessFunction(*points[i]) << '\n';
        delete[] points;
        points = nullptr;
	}
}