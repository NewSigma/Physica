/*
 * Copyright 2024 Weibo He.
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
#include <QApplication>
#include "Physica/Core/IO/Mnist.h"
#include "Physica/Core/AI/NeuralNetwork/Layer/LinearLayer.h"
#include "Physica/Core/AI/NeuralNetwork/SimpleNet.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Math/Optimization/Stochastic/MomentumSGD.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;
using namespace Physica::Gui;

template<Scalar> class MnistNet;

namespace Physica {
    template<class T>
    class Traits<MnistNet<T>> : public Traits<LinearLayer<T>> {};
}

template<Scalar T>
class MnistNet : public SimpleNet<MnistNet<T>> {
    using Base = SimpleNet<MnistNet<T>>;
    using typename Base::ValueType;
    using typename Base::InputType;
    using typename Base::OutputType;
private:
    LinearLayer<T> layer1;
    LinearLayer<T, false> layer2;
public:
    MnistNet() = default;
    template<RandomGenerator R>
    MnistNet(size_t width1, R&)
            : layer1(Mnist::NumPixelInImage, width1)
            , layer2(width1, 10) {
        //auto dist = std::normal_distribution<float>(0, 0.01);
        //layer1.template random_any<R>(dist);
        //layer2.template random_any<R>(dist);

        //layer1.template random_xavier_uniform<R>(1);
        //layer2.template random_xavier_uniform<R>(1);

        layer1.template random_xavier_normal<R>(1);
        layer2.template random_xavier_normal<R>(1);
    }
    template<Scalar U>
    MnistNet(const MnistNet<U>& net) : layer1(net.layer1), layer2(net.layer2) {}
    MnistNet(const MnistNet& other) = default;
    MnistNet(MnistNet&&) noexcept = default;
    ~MnistNet() = default;
    /* Operators */
    MnistNet& operator=(MnistNet& obj) noexcept { swap(obj); return *this; }
    /* Operations */
    [[nodiscard]] OutputType forward(const InputType& x) const {
        OutputType result = relu(layer1.forward(x));
        result = layer2.forward(result);
        return result;
    }

    template<class Dataset>
    [[nodiscard]] T loss(const Dataset& dataset, size_t index) const {
        return Loss<T>::crossEntropy(forward(dataset.getSamples()[index]), dataset.getLabels()[index]);
    }

    template<class Dataset>
    [[nodiscard]] T loss(const Dataset& dataset) const { return Base::loss(dataset); }

    template<class Dataset>
    ValueType calcAccuracy(const Dataset& dataset) const {
        const auto& testSamples = dataset.getSamples();
        const auto& testLabels = dataset.getLabels();
        const size_t numTestData = dataset.getSize();
        size_t count = 0;
        for (size_t i = 0; i < numTestData; ++i)
            count += testLabels[i] == Base::classify(InputType(testSamples[i]));
        return ValueType(count) / ValueType(numTestData);
    }

    void swap(MnistNet& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        layer1.swap(obj.layer1);
        layer2.swap(obj.layer2);
    }
private:
    template<Scalar> friend class MnistNet;
};

using ValueType = float32;
using ScalarType = Diff<ValueType, DiffMode::Reverse>;
using Dataset = Mnist::DatasetType<VectorND<ValueType>>;
using Optimizer = MomentumSGD<ScalarType>;
using RandomType = Random<MT19937>;
constexpr size_t numEpoch = 10;
constexpr size_t batchSize = 64;
constexpr double momentum = 0.5;
constexpr double learnRate = 0.05;

std::pair<Dataset, Dataset> makeDataset() {
    const Mnist mnist("/home/sigma/Documents/data");
    auto dataset = mnist.makeTrainDataset<VectorND<ValueType>>();
    for (size_t i = 0; i < dataset.getSize(); ++i) {
        auto& sample = dataset.getSamples()[i];
        sample = sample * ValueType(1.0 / 128) - ValueType(1);
    }
    return dataset.randomSplit<RandomType>(54000);
}

int main(int argc, char** argv) {
    ThreadPool::numThreadRequired = 4;
    const auto dataset = makeDataset();
    const size_t itePerEpoch = (dataset.first.getSize() + batchSize - 1) / batchSize;

    auto opt = Optimizer(momentum, learnRate, batchSize);
    opt.recordBegin();
    using NetType = MnistNet<ScalarType>;
    auto nn = NetType(32, RandomType::getInstance());
    opt.recordEnd();

    VectorND<ValueType> loss_train(numEpoch), loss_valid(numEpoch), acc_train(numEpoch), acc_valid(numEpoch);
    for (size_t epoch = 0; epoch < numEpoch; ++epoch) {
        if (epoch != 0) {
            for (size_t i = 0; i < itePerEpoch; ++i)
                nn.train_step<Dataset, Optimizer, RandomType, ThreadExecutor>(dataset.first, opt);
        }

        const auto nn_infer = MnistNet<ValueType>(nn);
        loss_train[epoch] = nn_infer.loss(dataset.first);
        loss_valid[epoch] = nn_infer.loss(dataset.second);
        acc_train[epoch] = nn_infer.calcAccuracy(dataset.first);
        acc_valid[epoch] = nn_infer.calcAccuracy(dataset.second);
    }

    QApplication app(argc, argv);
    Plot* plot = new Plot(-0.5, numEpoch - 0.5, 0, 3, 2, 1);
    auto* legend = plot->getChart()->legend();
    legend->setAlignment(Qt::AlignRight);
    legend->setMarkerShape(QLegend::MarkerShapeFromSeries);
    auto* axisX = plot->getAxisX();
    auto* axisY = plot->getAxisY();
    auto* axisRight = plot->getAxisRight();
    axisX->setTitleText("Epoch");
    axisX->setLabelFormat("%d");
    axisY->setTitleText("Loss");
    axisY->setLabelFormat("%d");
    axisRight->setTitleText("Accuracy");
    axisRight->setLabelsVisible(true);
    axisRight->setLabelFormat("%.1f");
    axisRight->setRange(0, 1);
    axisRight->setTickInterval(0.2);
    {
        auto& line = plot->line(loss_train);
        auto pen = line.pen();
        pen.setColor(Qt::black);
        line.setPen(pen);
        line.setName("Loss(train)");
    }
    {
        auto& line = plot->line(loss_valid);
        auto pen = line.pen();
        pen.setColor(Qt::red);
        line.setPen(pen);
        line.setName("Loss(valid)");
    }
    {
        auto& line = plot->line(acc_train);
        auto pen = line.pen();
        pen.setColor(Qt::black);
        pen.setStyle(Qt::DashLine);
        line.setPen(pen);
        line.setName("Acc(train)");
        line.detachAxis(axisY);
        line.attachAxis(axisRight);
    }
    {
        auto& line = plot->line(acc_valid);
        auto pen = line.pen();
        pen.setColor(Qt::red);
        pen.setStyle(Qt::DashLine);
        line.setPen(pen);
        line.setName("Acc(valid)");
        line.detachAxis(axisY);
        line.attachAxis(axisRight);
    }
    plot->show();
    return QApplication::exec();
}
