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
#include <Physica/Core/IO/Mnist.h>
#include <Physica/Core/AI/NeuralNetwork/Layer/LinearLayer.cuh>
#include <Physica/Core/AI/NeuralNetwork/SimpleNet.cuh>
#include <Physica/Core/AI/NeuralNetwork/SimpleDataset.h>
#include <Physica/Core/Math/Random/Random.h>
#include <Physica/Core/Math/Optimization/Stochastic/SGD.cuh>
#include <Physica/Gui/Plot/Plot.h>

using namespace Physica::Core;
using namespace Physica::Gui;

template<class ScalarType> class MnistNet;

namespace Physica {
    template<class T>
    class Traits<Core::device_obj<MnistNet<T>>> : public Traits<Core::device_obj<LinearLayer<T>>> {};
}

namespace Physica::Core {
    template<class ScalarType>
    class device_obj<MnistNet<ScalarType>> : public device_obj<SimpleNet<MnistNet<ScalarType>>> {
        using host_obj = MnistNet<ScalarType>;
        using This = device_obj<host_obj>;
        using Base = device_obj<SimpleNet<host_obj>>;
        using typename Base::ValueType;
        using typename Base::InputType;
        using typename Base::OutputType;
        using typename Base::LossType;
        using Base::IsTrainMode;
    public:
        using device_obj_type = This;
    private:
        device_obj<LinearLayer<ScalarType>> layer1;
        device_obj<LinearLayer<ScalarType, false>> layer2;
    public:
        device_obj() = default;
        template<class RandomGenerator>
        device_obj(size_t width1, RandomGenerator& gen)
                : layer1(decltype(layer1)::random_xavier_normal(Mnist::NumPixelInImage, width1, 1, gen))
                , layer2(decltype(layer2)::random_xavier_normal(width1, 10, 1, gen)) {}
        template<class OtherScalar>
        device_obj(const device_obj<MnistNet<OtherScalar>>& net) : layer1(net.getLayer1()), layer2(net.getLayer2()) {}
        device_obj(const This& other) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(This& obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] OutputType forward(const InputType& x) const {
            OutputType result = relu(layer1.forward(x));
            result = layer2.forward(result);
            return result;
        }

        template<class Dataset>
        [[nodiscard]] LossType loss(const Dataset& dataset, size_t index) const {
            OutputType output;
            if constexpr (IsTrainMode)
                output = forward(dataset.getSamples()[index]);
            else
                output = forward(dataset.getSamples()[index].getValues());
            return device_obj<Loss<ScalarType>>::crossEntropy(output, dataset.getLabels()[index]);
        }

        template<class Dataset>
        [[nodiscard]] LossType loss(const Dataset& dataset) const { return Base::loss(dataset); }

        [[nodiscard]] This copy() const {
            This result{};
            result.layer1 = layer1.copy();
            result.layer2 = layer2.copy();
            return result;
        }

        template<class Dataset>
        ValueType calcAccuracy(const Dataset& dataset) const {
            const auto& testSamples = dataset.getSamples();
            const auto& testLabels = dataset.getLabels();
            const size_t numTestData = dataset.getSize();
            size_t count = 0;
            for (size_t i = 0; i < numTestData; ++i)
                count += testLabels[i] == Base::classify(testSamples[i]);
            return ValueType(count) / ValueType(numTestData);
        }

        void swap(This& __restrict obj) noexcept {
            assert(this != &obj && "[Error]: Self swap is likely a bug");
            Base::swap(obj);
            layer1.swap(obj.layer1);
            layer2.swap(obj.layer2);
        }
        /* Getters */
        [[nodiscard]] const auto& getLayer1() const noexcept { return layer1; }
        [[nodiscard]] const auto& getLayer2() const noexcept { return layer2; }
    };
}

using ValueType = float32;
using ScalarType = Diff<ValueType, DiffMode::Reverse, 1>;
using DeviceVector = device_obj<Diff<VectorND<ValueType>, DiffMode::Reverse, 1>>;
using Dataset = typename Mnist::DatasetType<DeviceVector>;
using Optimizer = SGD<device_obj<ScalarType>>;
using RandomGenerator = std::mt19937;
using RandomType = Random<RandomGenerator>;
constexpr size_t numEpoch = 10;
constexpr size_t batchSize = 64;
constexpr double learnRate = 0.05;

std::pair<Dataset, Dataset> makeDataset(RandomGenerator& gen) {
    const Mnist mnist("/home/sigma/Documents/data");
    auto dataset = mnist.makeTrainDataset<VectorND<ValueType>>();
    for (size_t i = 0; i < dataset.getSize(); ++i) {
        auto& sample = dataset.getSamples()[i];
        sample = sample * ValueType(1.0 / 128) - ValueType(1);
    }
    return dataset.randomSplit(54000, gen);
}

int main(int argc, char** argv) {
    ThreadPool::numThreadRequired = 4;
    auto& gen = RandomType::getInstance().getGen();
    const auto dataset = makeDataset(gen);
    const size_t itePerEpoch = (dataset.first.getSize() + batchSize - 1) / batchSize;

    auto opt = Optimizer(learnRate, batchSize);
    opt.recordBegin();
    auto nn = device_obj<MnistNet<ScalarType>>(32, gen);
    opt.recordEnd();

    VectorND<ValueType> loss_train(numEpoch), loss_valid(numEpoch), acc_train(numEpoch), acc_valid(numEpoch);
    for (size_t epoch = 0; epoch < numEpoch; ++epoch) {
        if (epoch != 0) {
            for (size_t i = 0; i < itePerEpoch; ++i)
                nn.train_step<Dataset, Optimizer, RandomType>(dataset.first, opt);
        }

        using VectorType = VectorND<ValueType>;
        const auto nn_infer = device_obj<MnistNet<ValueType>>(nn);
        loss_train[epoch] = nn_infer.loss(dataset.first);
        loss_valid[epoch] = nn_infer.loss(dataset.second);
        acc_train[epoch] = nn_infer.calcAccuracy(dataset.first);
        acc_valid[epoch] = nn_infer.calcAccuracy(dataset.second);
    }

    QApplication app(argc, argv);
    Plot* plot = new Plot(-0.5, numEpoch - 0.5, 0, 3, 2, 1);
    auto& legend = plot->getLegend();
    legend.setAlignment(Qt::AlignRight);
    legend.show();
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
