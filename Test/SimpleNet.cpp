/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <Core/Header/Vector.h>
#include <AI/Header/DNN.h>
#include <AI/Header/Layer.h>
#include <AI/Header/Node.h>

void simpleNet() {
    int* arr = new int[3]{2, 2, 2};
    auto& net = *new Physica::AI::DNN(2, 3, arr);
    net.setLearnRate(Numerical(0.1));
    net[0][0].connect(0, 0);
    net[0][0].connect(0, 1);
    net[0][1].connect(1, 0);
    net[0][1].connect(1, 2);
    net[1][0].connect(0, 0);
    net[1][0].connect(0, 2);
    net[1][1].connect(1, 0);
    net[1][1].connect(1, 1);
    Vector d1(Numerical((SignedNumericalUnit)1), Numerical((SignedNumericalUnit)1));
    Numerical n1((SignedNumericalUnit)1);
    net.loadData(d1, n1);
    for(int i = 0; i < 5; ++i) {
        net.train();
        std::cout << "Train " << i << " finished. Loss:" << double(net.predict()) << '\n';
    }
    Vector d2(Numerical((SignedNumericalUnit)1.1), Numerical((SignedNumericalUnit)0.9));
    Numerical n2((SignedNumericalUnit)1);
    net.loadData(d2, n2);
    std::cout << "Test. Loss:" << double(net.predict()) << '\n';
    delete[] arr;
    delete &net;
}