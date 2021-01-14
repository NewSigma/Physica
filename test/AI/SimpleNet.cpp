/*
 * Copyright 2019 WeiBo He.
 *
 * This file is part of Physica.
 *
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
#include <iostream>
#include <Physica/Physica.h>
#include <Physica/AI/DNN.h>
#include <Physica/AI/Layer.h>
#include <Physica/AI/Node.h>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector.h>

using namespace Physica::Core;

int main(int argc, char** argv) {
    initPhysica();
    int* arr = new int[3]{2, 2, 2};
    auto& net = *new Physica::AI::DNN(2, 3, arr);
    net.setLearnRate(MultiScalar(0.1));
    net[0][0].connect(0, 0);
    net[0][0].connect(0, 1);
    net[0][1].connect(1, 0);
    net[0][1].connect(1, 2);
    net[1][0].connect(0, 0);
    net[1][0].connect(0, 2);
    net[1][1].connect(1, 0);
    net[1][1].connect(1, 1);
    Vector<MultiScalar> d1(CStyleArray<MultiScalar>(2));
    d1 << 1.0 << 1.0;
    MultiScalar n1((SignedScalarUnit)1);
    net.loadData(d1, n1);
    for(int i = 0; i < 5; ++i) {
        net.train();
        std::cout << "Train " << i << " finished. Loss:" << double(net.predict()) << '\n';
    }
    Vector<MultiScalar> d2(CStyleArray<MultiScalar>(2));
    d2 << 1.1 << 0.9;
    MultiScalar n2((SignedScalarUnit)1);
    net.loadData(d2, n2);
    std::cout << "test. Loss:" << double(net.predict()) << '\n';
    delete[] arr;
    delete &net;
    deInitPhysica();
    return 0;
}