/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector.h>
#include <Physica/AI/DNN.h>
#include <Physica/AI/Layer.h>
#include <Physica/AI/Node.h>
#include <iostream>

namespace Physica::Test {
    using namespace Physica::Core;
    void simpleNet() {
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
        Vector<MultiScalar> d1(CStyleArray<MultiScalar, Dynamic, Dynamic>(2));
        d1 << 1 << 1;
        MultiScalar n1((SignedScalarUnit)1);
        net.loadData(d1, n1);
        for(int i = 0; i < 5; ++i) {
            net.train();
            std::cout << "Train " << i << " finished. Loss:" << double(net.predict()) << '\n';
        }
        Vector<MultiScalar> d2(CStyleArray<MultiScalar, Dynamic, Dynamic>(2));
        d2 << 1.1 << 0.9;
        MultiScalar n2((SignedScalarUnit)1);
        net.loadData(d2, n2);
        std::cout << "test. Loss:" << double(net.predict()) << '\n';
        delete[] arr;
        delete &net;
    }
}