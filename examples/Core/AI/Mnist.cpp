/*
 * Copyright 2024 WeiBo He.
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
#include "Physica/Core/AI/NeuralNetwork/NetBase.h"
#include "Physica/Core/AI/NeuralNetwork/Loss.h"
#include "Physica/Core/Math/Random/RandomPool.h"
#include "Physica/Core/Math/Optimization/Stochastic/MomentumSGD.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;
using namespace Physica::Gui;

template<class ScalarType> class MnistNet;

namespace Physica::Core::Internal {
    template<class T>
    class Traits<MnistNet<T>> : public Traits<LinearLayer<T>> {};
}

template<class ScalarType>
class MnistNet : public NetBase<MnistNet<ScalarType>> {
    using Base = NetBase<MnistNet<ScalarType>>;
    using typename Base::PlainScalar;
    using typename Base::InputType;
    using typename Base::OutputType;
private:
    LinearLayer<ScalarType> layer1;
    LinearLayer<ScalarType, false> layer2;
public:
    MnistNet() = default;
    template<class RandomGenerator>
    MnistNet(size_t width1, RandomGenerator& gen)
            : layer1(Mnist::NumPixelInImage, width1)
            , layer2(width1, 10) {
        //auto dist = std::normal_distribution<float>(0, 0.01);
        //layer1.random_any(dist, gen);
        //layer2.random_any(dist, gen);

        //layer1.random_xavier_uniform(1, gen);
        //layer2.random_xavier_uniform(1, gen);

        layer1.random_xavier_normal(1, gen);
        layer2.random_xavier_normal(1, gen);
    }
    template<class OtherScalar>
    MnistNet(const MnistNet<OtherScalar>& net) : layer1(net.layer1), layer2(net.layer2) {}
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
    [[nodiscard]] ScalarType loss(const Dataset& dataset, size_t index) const {
        return Loss<ScalarType>::crossEntropy(forward(dataset.getSamples()[index]), dataset.getLabels()[index]);
    }

    template<class Dataset>
    [[nodiscard]] ScalarType loss(const Dataset& dataset) const { return Base::loss(dataset); }

    [[nodiscard]] MnistNet copy() const {
        MnistNet result{};
        result.layer1 = layer1.copy();
        result.layer2 = layer2.copy();
        return result;
    }

    template<class Dataset>
    PlainScalar calcAccuracy(const Dataset& dataset) const {
        const auto& testSamples = dataset.getSamples();
        const auto& testLabels = dataset.getLabels();
        const size_t numTestData = dataset.getSize();
        size_t count = 0;
        for (size_t i = 0; i < numTestData; ++i)
            count += testLabels[i] == Base::classify(InputType(testSamples[i]));
        return PlainScalar(count) / PlainScalar(numTestData);
    }

    void swap(MnistNet& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        layer1.swap(obj.layer1);
        layer2.swap(obj.layer2);
    }
private:
    template<class OtherScalar> friend class MnistNet;
};

using PlainScalar = Scalar<Float>;
using ScalarType = Differentiable<PlainScalar, DiffMode::Reverse>;
using Dataset = typename Mnist::DatasetType<PlainScalar>;
using Optimizer = MomentumSGD<ScalarType>;
using RandomGenerator = std::mt19937;
using RandomPoolType = RandomPool<RandomGenerator>;
constexpr size_t numEpoch = 50;
constexpr size_t batchSize = 64;
constexpr double momentum = 0.5;
constexpr double learnRate = 0.05;

std::pair<Dataset, Dataset> makeDataset(RandomGenerator& gen) {
    const Mnist mnist("/home/sigma/Documents/data");
    auto dataset = mnist.makeTrainDataset<PlainScalar>();
    for (size_t i = 0; i < dataset.getSize(); ++i) {
        auto& sample = dataset.getSamples()[i];
        sample = sample * PlainScalar(1.0 / 128) - PlainScalar(1);
    }
    return dataset.randomSplit(54000, gen);
}

int main(int argc, char** argv) {
    ThreadPool::numThreadRequired = 4;
    auto& gen = RandomPoolType::getInstance().getGen();
    const auto dataset = makeDataset(gen);
    const size_t itePerEpoch = (dataset.first.getSize() + batchSize - 1) / batchSize;

    auto opt = Optimizer(momentum, learnRate, batchSize);
    opt.recordBegin();
    auto nn = MnistNet<ScalarType>(32, gen);
    opt.recordEnd();

    Vector<PlainScalar> loss_train(numEpoch), loss_valid(numEpoch), err_train(numEpoch), err_valid(numEpoch);
    for (size_t epoch = 0; epoch < numEpoch; ++epoch) {
        for (size_t i = 0; i < itePerEpoch; ++i)
            nn.train_step<Dataset, Optimizer, RandomPoolType, ThreadExecutor>(dataset.first, opt);
        const auto nn_infer = MnistNet<PlainScalar>(nn);
        loss_train[epoch] = nn_infer.loss(dataset.first);
        loss_valid[epoch] = nn_infer.loss(dataset.second);
        err_train[epoch] = PlainScalar(1) - nn_infer.calcAccuracy(dataset.first);
        err_valid[epoch] = PlainScalar(1) - nn_infer.calcAccuracy(dataset.second);
    }

    QApplication app(argc, argv);
    Plot* plot = new Plot(-0.5, numEpoch - 0.5, 0, 0.35, 20, 0.1);
    auto* legend = plot->getChart()->legend();
    legend->setAlignment(Qt::AlignTop);
    legend->setMarkerShape(QLegend::MarkerShapeFromSeries);
    auto* axisX = plot->getAxisX();
    auto* axisY = plot->getAxisY();
    auto* axisRight = plot->getAxisRight();
    axisX->setTitleText("Epoch");
    axisX->setLabelFormat("%d");
    axisY->setTitleText("Loss");
    axisY->setLabelFormat("%.1f");
    axisRight->setTitleText("Error rate");
    axisRight->setLabelsVisible(true);
    axisRight->setLabelFormat("%.1f");
    {
        auto& line = plot->line(loss_train);
        auto pen = line.pen();
        pen.setColor(Qt::black);
        line.setPen(pen);
        line.setName("loss_train");
    }
    {
        auto& line = plot->line(loss_valid);
        auto pen = line.pen();
        pen.setColor(Qt::red);
        line.setPen(pen);
        line.setName("loss_valid");
    }
    {
        auto& line = plot->line(err_train);
        auto pen = line.pen();
        pen.setColor(Qt::black);
        pen.setStyle(Qt::DashLine);
        line.setPen(pen);
        line.setName("err_train");
        line.detachAxis(axisY);
        line.attachAxis(axisRight);
    }
    {
        auto& line = plot->line(err_valid);
        auto pen = line.pen();
        pen.setColor(Qt::red);
        pen.setStyle(Qt::DashLine);
        line.setPen(pen);
        line.setName("err_valid");
        line.detachAxis(axisY);
        line.attachAxis(axisRight);
    }
    plot->show();
    return QApplication::exec();
}
