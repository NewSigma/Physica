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
#ifndef PHYSICA_TOTALCIRCUIT_H
#define PHYSICA_TOTALCIRCUIT_H

#include <vector>
#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Utils/Container/FreezingList/FreezingList.h"

namespace Physica::Core {
    class TotalCircuit {
        using ScalarType = Scalar<Double, false>;
        /*
         * Save two connection points.
         *
         * Note: We always assume pos1 < pos2 for the convenience of implementation.
         */
        struct Connection {
        private:
            ScalarType* resistance;
            int pos1, pos2;
        public:
            inline void setData(ScalarType* p_scalar, int i1, int i2) {
                assert(i1 < i2);
                resistance = p_scalar;
                pos1 = i1;
                pos2 = i2;
            }
            [[nodiscard]] inline ScalarType* getResistance() { return resistance; }
            [[nodiscard]] inline int getPos1() const { return pos1; }
            [[nodiscard]] inline int getPos2() const { return pos2; }
        };
    private:
        //Value of resistances.
        Utils::FreezingList<ScalarType> resistances;
        std::vector<ScalarType> result;
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
        explicit TotalCircuit(std::vector<ScalarType> v);
        TotalCircuit(const TotalCircuit&) = delete;
        TotalCircuit(TotalCircuit&&) noexcept = delete;
        ~TotalCircuit();
        /* Operators */
        TotalCircuit& operator=(const TotalCircuit&) = delete;
        TotalCircuit& operator=(TotalCircuit&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] const std::vector<ScalarType>& getResult() const noexcept { return result; }
    private:
        void connect();
        void calculate();
    };
}


#endif