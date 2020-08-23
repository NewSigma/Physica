/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef _Physica_C_GeneAlgorithm_H
#define _Physica_C_GeneAlgorithm_H

#include <ctime>

#include "Physica/Core/MultiPrecition/Scalar.h"

namespace Physica::Core {
    class GeneAlgorithm {
    public:
        enum ChooseMode {
            LinearChoose,
            RandomChoose
        };
        double crossoverRate = 0.6;
        double mutationRate = 0.1;

        GeneAlgorithm(MultiScalar func(const MultiScalar&), const MultiScalar& lower, const MultiScalar& upper, int pop = 100, ChooseMode mode = RandomChoose);
        ~GeneAlgorithm();

        MultiScalar** getExtremalPoint();
        void setMaxGenerations(int maxGenerations);
        void setMaxTime(int maxTime);
        void print();
    private:
        //Algorithm basic args.
        //The size of initial points we should choose.
        int population;
        //Function args.
        MultiScalar (*fitnessFunction)(const MultiScalar&);
        const MultiScalar* lowerBound;
        const MultiScalar* upperBound;
        MultiScalar* regionLength;
        //Save
        MultiScalar** points;
        //Stopping args
        int generations = 0;
        clock_t startTime{};
        int maxGenerations = 100;
        int maxTime = 600000;

        void crossover();
        void mutation();
        bool shouldStop();
    };
}



#endif