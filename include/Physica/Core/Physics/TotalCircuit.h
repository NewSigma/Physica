/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_TOTALCIRCUIT_H
#define PHYSICA_TOTALCIRCUIT_H

#include <vector>
#include "Physica/Core/MultiPrecition/Scalar.h"
#include "Physica/Core/Utils/FreezingList/FreezingList.h"

namespace Physica::Core {
    class TotalCircuit {
        /*
         * Save two connection points.
         *
         * Note: We always assume pos1 < pos2 for the convenience of implementation.
         */
        struct Connection {
        private:
            MultiScalar* resistance;
            int pos1, pos2;
        public:
            inline void setData(MultiScalar* p_scalar, int i1, int i2) {
                Q_ASSERT(i1 < i2);
                resistance = p_scalar;
                pos1 = i1;
                pos2 = i2;
            }
            [[nodiscard]] inline MultiScalar* getResistance() { return resistance; }
            [[nodiscard]] inline int getPos1() const { return pos1; }
            [[nodiscard]] inline int getPos2() const { return pos2; }
        };
    private:
        //Value of resistances.
        FreezingList<MultiScalar> resistances;
        std::vector<MultiScalar> result;
        //Size of all resistances. Not always equals to resistances.size().
        size_t size;
        //Total number of nodes.
        int nodesCount;
        /*
         * Temp array, store order produced by solve(),
         * we may pass it to calculate the equivalent resistance or simply save it.
         */
        Connection* order;
        //Save the first element of order that is not initialized.
        Connection* orderCurrent;
    public:
        explicit TotalCircuit(std::vector<MultiScalar> v);
        TotalCircuit(const TotalCircuit&) = default;
        TotalCircuit(TotalCircuit&&) noexcept = default;
        ~TotalCircuit();
        /* Operators */
        TotalCircuit& operator=(const TotalCircuit&) = default;
        TotalCircuit& operator=(TotalCircuit&&) noexcept = default;
        /* Getters */
        [[nodiscard]] const std::vector<MultiScalar>& getResult() const noexcept { return result; }
    private:
        void connect();
        void calculate();
    };
}


#endif