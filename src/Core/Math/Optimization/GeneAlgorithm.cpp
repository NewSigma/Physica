/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <iostream>
#include "Physica/Core/Math/Optimization/GeneAlgorithm.h"
#include "Physica/Core/ElementaryFunction.h"
#include "Physica/Core/Scalar.h"

namespace Physica::Core {
    //Usage: GeneAlgorithm(args)->initFunction(args)->getExtremalPoint()
    GeneAlgorithm::GeneAlgorithm(Scalar func(const Scalar&), const Scalar& lower, const Scalar& upper, int pop, ChooseMode mode) {
        population = pop;
        points = new Scalar*[population];

        fitnessFunction = func;
        lowerBound = new Scalar(lower);
        upperBound = new Scalar(upper);
        regionLength = new Scalar(upper - lower);
        //Get abs(regionLength).
        regionLength->toAbs();

        if (mode == LinearChoose) {
            Scalar stepLength = *regionLength / Scalar(static_cast<SignedScalarUnit>(population));
            Scalar temp(1, 0);
            for (int i = 0; i < population; i++) {
                temp[0] = i;
                points[i] = new Scalar(stepLength * temp + lower);
            }
        }
        else if (mode == RandomChoose) {
            for (int i = 0; i < population; i++)
                points[i] = new Scalar(randomNumerical() * *regionLength + lower);
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
    Scalar** GeneAlgorithm::getExtremalPoint() {
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
            if (Scalar(crossoverRate) > r) {
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
        if (Scalar(mutationRate) > randomNumerical()) {
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
}