#include "GeneAlgorithm.h"
#include <iostream>
#include "ElementaryFunction.h"
#include "Numerical.h"
/*
 * Usage: GeneAlgorithm(args)->initFunction(args)->getExtremalPoint()
 * Copyright (c) 2019 NewSigma@163.com.All rights reserved.
 */
//TODO Debug
GeneAlgorithm::GeneAlgorithm(Numerical* func(const Numerical&), const Numerical* lower, const Numerical* upper, int pop, ChooseMode mode) {
    population = pop;
	points = new Numerical*[population];

    fitnessFunction = func;
    lowerBound = lower;
    upperBound = upper;
    regionLength = *upperBound - *lowerBound;
    //Get abs(regionLength).
    regionLength->length = (signed char)regionLength->getSize();

    Numerical* element;
    if (mode == LinearChoose) {
        auto stepLength = *regionLength / Numerical(population);
        Numerical* number_i;
        for (int i = 0; i < population; i++) {
            number_i = new Numerical(i);
            element = *stepLength * *number_i;
            *element += *lowerBound;
            points[i] = element;
            delete number_i;
        }
        delete stepLength;
    }
    else if (mode == RandomChoose) {
        for (int i = 0; i < population; i++) {
            element = new Numerical(randomNumerical());
            *element *= *regionLength;
            *element += *lowerBound;
            points[i] = element;
        }
    }
}

GeneAlgorithm::~GeneAlgorithm() {
    delete[] points;
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
		if (Numerical(crossoverRate) > *r) {
			long randomIndex1 = random() % population;
			long randomIndex2 = random() % population;
            auto random1 = points[randomIndex1];
            auto random2 = points[randomIndex2];

            delete r;
			r = randomNumerical();
            //Whether random2 > random1 or not is not important.
            auto child = *random2 - *random1;
            *child *= *r;
            *child += *random1;
            auto child_y = fitnessFunction(*child);
            auto y_random1 = fitnessFunction(*random1);
            auto y_random2 = fitnessFunction(*random2);
			if (*child_y > *y_random1)
                *random1 = *child;
			else if (*child_y > *y_random2)
                *random2 = *child;

			delete child;
			delete child_y;
			delete y_random1;
			delete y_random2;
		}
		delete r;
	}
}

void GeneAlgorithm::mutation() {
    Numerical* r = randomNumerical();
	if (Numerical(mutationRate) > *r) {
		long randomIndex = random() % population;
		*points[randomIndex] = *randomNumerical();
        *points[randomIndex] *= *regionLength;
        *points[randomIndex] += *lowerBound;
	}
	delete r;
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
        for (int i = 0; i < population; ++i) {
            Numerical* value = fitnessFunction(*points[i]);
            std::cout << *points[i] << " " << *value << '\n';
            delete value;
        }
        delete[] points;
        points = nullptr;
	}
}