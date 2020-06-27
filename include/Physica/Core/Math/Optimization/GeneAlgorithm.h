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