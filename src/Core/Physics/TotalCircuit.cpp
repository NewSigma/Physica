/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "Physica/Core/Physics/TotalCircuit.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/LinearEquations.h"

namespace Physica::Core {
    /*!
     * @param v is the list that contains the value of resistances.
     *
     * Initially nodesCount = 2, that is,
     * two terminals of the two terminals network, whose order numbers are 0 and (size + 1).
     *
     * Improve:
     * 1. The resistance must be positive, may be we should add a pack class to avoid incorrect use.
     * 2. Some results may be very large, which means the circuit is cut off, we do not interested in them,
     * try to remove them.
     * 3. Some results may be very same, we only have to save one of them.
     *
     * Optimize:
     * If two resistances are equal, they are changeable and the program will do useless work.
     */
    TotalCircuit::TotalCircuit(std::vector<MultiScalar> v) : resistances(v.size()), nodesCount(2) {
        //In case overflow
        size = resistances.size();
        Q_ASSERT(size < static_cast<size_t>(INT_MAX));

        for(auto& s : v)
            resistances.allocate(std::move(s));

        orderCurrent = order = new Connection[size];
        //Start
        connect();
    }

    TotalCircuit::~TotalCircuit() {
        delete[] order;
    }
    /*!
     * Using DFS, recursively select resistances.
     * When all the resistances is marked used, calculate the equivalent resistance.
     */
    void TotalCircuit::connect() {
        if(resistances.empty()) {
            calculate();
            return;
        }

        typedef FreezingList<MultiScalar>::Node Node;
        const auto end = resistances.end();
        for(auto ite = resistances.begin(); ite != end; ++ite) {
            resistances.freeze(ite);
            /*
             * We first handle the connection between node[0] and the other nodes, except empty node.
             * In the mean time, we will get the position of empty node.
             */
            const int nodesCount_1 = nodesCount - 1;
            Connection* copyCurrent = orderCurrent;
            ++orderCurrent;
            /* Handle node[0] */ {
                const int size_1 = static_cast<int>(size) + 1;
                for(int i = 0; i < nodesCount_1; ++i) {
                    copyCurrent->setData(&*ite, i, size_1);
                    connect();
                }
            }
            /* Handle other nodes */ {
                const int nodesCount_2 = nodesCount_1 - 1;
                for(int i = 0; i < nodesCount_2; ++i) {
                    for(int j = i + 1; j < nodesCount_1; ++j) {
                        copyCurrent->setData(&*ite, i, j);
                        connect();
                    }
                }
            }
            //Finally, we handle the empty node which connect to other nodes.
            ++nodesCount;
            for(int i = 0; i < nodesCount_1; ++i) {
                copyCurrent->setData(&*ite, i, nodesCount_1);
                connect();
            }
            --nodesCount;

            --orderCurrent;
            //This branch finished. Add the temp back.
            resistances.restore(ite);
        }
    }
    /*!
     * Using Kirchoff equations to calculate the equivalent resistance.
     */
    void TotalCircuit::calculate() {
        /*
         * Assume the voltage at node[0] is 0 and voltage at node[1] is 1.
         * We will have to solve a (nodesCount - 2) rank linear Equations to get the other voltages.
         */
        const int rank = nodesCount - 2;
        if(rank) {
            Matrix<MultiScalar, Row> augmentedMatrix;
            /* Construct augmented matrix */ {
                augmentedMatrix = Matrix<MultiScalar, Row>(rank);
                for(int i = 0; i < rank; ++i)
                    augmentedMatrix.allocate(Matrix<MultiScalar, Row>::VectorType::zeroVector(rank + 1), i);

                Connection* p = order;
                const int size_1 = static_cast<int>(size) + 1;
                while(p != orderCurrent) {
                    const int pos1 = p->getPos1();
                    const int pos2 = p->getPos2();

                    auto conductance = reciprocal(*p->getResistance());

                    if(pos2 != size_1) {
                        //index = pos - 1. Here the possible pos starts from 1.
                        const int index2 = pos2 - 1;
                        auto& row2 = augmentedMatrix[index2];
                        if(pos1 != 0) {
                            const int index1 = pos1 - 1;
                            auto& row1 = augmentedMatrix[index1];
                            row1[index1] += conductance;
                            row1[index2] -= conductance;
                            row2[index1] -= conductance;
                            row2[index2] += conductance;
                        }
                        else {
                            row2[index2] += conductance;
                            row2[rank] += conductance;
                        }
                    }
                    else {
                        if(pos1 != 0) {
                            const int index1 = pos1 - 1;
                            auto& row1 = augmentedMatrix[index1];
                            row1[index1] += conductance;
                        }
                    }
                    ++p;
                }
            }
            /* Solve the equations */
            LinearEquations<MultiScalar, Row> le(std::move(augmentedMatrix));
            le.solve(AbstractLinearEquations::GaussEliminationPartial);
            /* Calculate equivalent resistance */
            Connection* p = order;
            MultiScalar totalCurrent = MultiScalar::getZero();
            while(p != orderCurrent) {
                //Calculate the current flow out of node 0.(Not node[0])
                if(p->getPos1() == 0) {
                    auto index2 = p->getPos2() - 1;
                    auto current = reciprocal(*p->getResistance());
                    if(index2 != size)
                        current -= current * le.getResult(index2);
                    totalCurrent += current;
                }
                ++p;
            }
            /*
             * Voltage between node[0] and node[1] is assumed to be 1.
             * Calculate reciprocal directly so we have the equivalent resistance.
             */
            if(!totalCurrent.isZero())
                result.push_back(reciprocal(totalCurrent));
        }
        else {
            Connection* p = order;
            MultiScalar totalCurrent = MultiScalar::getZero();
            while(p != orderCurrent) {
                auto current = reciprocal(*p->getResistance());
                totalCurrent += current;
                ++p;
            }
            result.push_back(reciprocal(totalCurrent));
        }
    }
}