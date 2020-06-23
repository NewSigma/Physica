/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <Core/Header/Vector.h>
#include <AI/Header/DNN.h>
#include <AI/Header/Layer.h>
#include <AI/Header/Node.h>

namespace Physica::Test {
    void simpleNet() {
        int* arr = new int[3]{2, 2, 2};
        auto& net = *new Physica::AI::DNN(2, 3, arr);
        net.setLearnRate(Scalar(0.1));
        net[0][0].connect(0, 0);
        net[0][0].connect(0, 1);
        net[0][1].connect(1, 0);
        net[0][1].connect(1, 2);
        net[1][0].connect(0, 0);
        net[1][0].connect(0, 2);
        net[1][1].connect(1, 0);
        net[1][1].connect(1, 1);
        Vector d1(Scalar((SignedScalarUnit)1), Scalar((SignedScalarUnit)1));
        Scalar n1((SignedScalarUnit)1);
        net.loadData(d1, n1);
        for(int i = 0; i < 5; ++i) {
            net.train();
            std::cout << "Train " << i << " finished. Loss:" << double(net.predict()) << '\n';
        }
        Vector d2(Scalar((SignedScalarUnit)1.1), Scalar((SignedScalarUnit)0.9));
        Scalar n2((SignedScalarUnit)1);
        net.loadData(d2, n2);
        std::cout << "test. Loss:" << double(net.predict()) << '\n';
        delete[] arr;
        delete &net;
    }
}