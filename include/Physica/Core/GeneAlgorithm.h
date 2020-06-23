#ifndef _Physica_C_GeneAlgorithm_H
#define _Physica_C_GeneAlgorithm_H

#include <ctime>

namespace Physica::Core {
    class Scalar;

    class GeneAlgorithm {
    public:
        enum ChooseMode {
            LinearChoose,
            RandomChoose
        };
        double crossoverRate = 0.6;
        double mutationRate = 0.1;

        GeneAlgorithm(Scalar func(const Scalar&), const Scalar& lower, const Scalar& upper, int pop = 100, ChooseMode mode = RandomChoose);
        ~GeneAlgorithm();

        Scalar** getExtremalPoint();
        void setMaxGenerations(int maxGenerations);
        void setMaxTime(int maxTime);
        void print();
    private:
        //Algorithm basic args.
        //The size of initial points we should choose.
        int population;
        //Function args.
        Scalar (*fitnessFunction)(const Scalar&);
        const Scalar* lowerBound;
        const Scalar* upperBound;
        Scalar* regionLength;
        //Save
        Scalar** points;
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