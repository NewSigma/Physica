/*
 * Copyright 2024-2025 Weibo He.
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
#include <QApplication>
#include "Physica/Core/IO/Mnist.h"
#include "Physica/Core/ML/NeuralNetwork/Layer/LinearLayer.h"
#include "Physica/Core/ML/NeuralNetwork/SeqNet.h"
#include "Physica/Core/ML/Optimizer/SGD.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Scalar/Diff.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica;
constexpr size_t NumEpoch = 10;
constexpr int BatchSize = 64;
constexpr int HiddenW = 32;
constexpr double LearnRate = 0.05;

namespace Physica {
    template<Scalar> class MnistNet;

    template<Scalar T>
    class Traits<MnistNet<T>> : public Traits<LinearLayer<T>> {};

    template<Scalar T>
    class MnistNet : public SeqNet<MnistNet<T>> {
        using This = MnistNet<T>;
        using Base = SeqNet<This>;
        using typename Base::Tv;
    private:
        LinearLayer<T> layer1 = LinearLayer<T>(Mnist::NumPixelInImage, HiddenW);
        LinearLayer<T> layer2 = LinearLayer<T>(HiddenW, 10);
        SGD<Tv> opt = SGD<Tv>(LearnRate);
    public:
        template<RNG R = Random<>>
        MnistNet(R&) {
            layer1.template random_xavier_normal<R>(1);
            layer2.template random_xavier_normal<R>(1);
        }
        template<Scalar U>
        MnistNet(const MnistNet<U>& net) : layer1(net.layer1), layer2(net.layer2) {}
        MnistNet(const This& other) = default;
        MnistNet(This&&) noexcept = default;
        ~MnistNet() = default;
        /* Operators */
        This& operator=(This& obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] CoDiff<VectorND<T>> forward(const Vector auto& x) const {
            auto y1 = layer1.forward(x);
            CoDiff<VectorND<T>> y2 = relu(y1);
            auto y3 = layer2.forward(y2);
            if constexpr (Base::IsTrain)
                y3 = co_yield std::move(y3);
            else
                co_return std::move(y3);
        }

        void reverse(const This& __restrict other) const noexcept {
            assert(this != &other && "[Error]: Self reverse is invalid");
            layer1.reverse(other.layer1);
            layer2.reverse(other.layer2);
        }

        void step(auto& optimizer) {
            layer1.step(optimizer);
            layer2.step(optimizer);
        }

        void step() { step(opt); }

        void zero_grad() {
            layer1.zero_grad();
            layer2.zero_grad();
        }

        [[nodiscard]] CoDiff<T> loss(const auto& dataset, size_t index) const {
            auto weights = forward(dataset.getSamples()[index]);
            auto loss = weights.crossEntropy(dataset.getLabels()[index]);
            if constexpr (ReverseDiff<T>) {
                auto tmp = co_yield loss.value();
                loss.reverse(tmp.grad());
            }
            else
                co_return std::move(loss);
        }

        [[nodiscard]] T loss(const auto& dataset) const { return Base::loss(dataset); }

        size_t classify(const VectorND<Tv>& input) const {
            static_assert(!Base::IsTrain, "[Error]: It is suggested using eval mode to reduce memory use");
            const auto output = forward(input);
            Tv max = output[0];
            size_t index = 0;
            for (size_t i = 1; i < output.getLength(); ++i) {
                if (output[i] > max) {
                    index = i;
                    max = output[i].value();
                }
            }
            return index;
        }

        Tv calcAccuracy(const auto& dataset) const {
            const auto& testSamples = dataset.getSamples();
            const auto& testLabels = dataset.getLabels();
            const size_t numTestData = dataset.getSize();
            size_t count = 0;
            for (size_t i = 0; i < numTestData; ++i)
                count += testLabels[i] == classify(VectorND<Tv>(testSamples[i]));
            return Tv(count) / Tv(numTestData);
        }

        void swap(This& __restrict obj) noexcept {
            assert(this != &obj && "[Error]: Self swap is likely a bug");
            Base::swap(obj);
            layer1.swap(obj.layer1);
            layer2.swap(obj.layer2);
        }
    private:
        template<Scalar> friend class MnistNet;
    };
}

using T = float32;
using dfloat = Diff<T, DiffMode::Reverse>;
using Dataset = Mnist::DatasetType<VectorND<T>>;
using RandomSource = Random<>;

namespace {
    std::pair<Dataset, Dataset> makeDataset() {
        const Mnist mnist("./data");
        auto dataset = mnist.makeTrainDataset<VectorND<T>>();
        for (size_t i = 0; i < dataset.getSize(); ++i) {
            auto& sample = dataset.getSamples()[i];
            sample = sample * T(1.0 / 128) - T(1);
        }
        return dataset.randomSplit<RandomSource>(54000);
    }
}

int main(int argc, char** argv) {
    ThreadPool::numThreadRequired = 4;

    const auto dataset = makeDataset();
    const int64_t stepPerEpoch = int64_t(dataset.first.getSize() + BatchSize - 1) / BatchSize;
    auto nn = MnistNet<dfloat>(RandomSource::getInstance());
    VectorND<T> loss_train(NumEpoch);
    VectorND<T> loss_valid(NumEpoch);
    VectorND<T> acc_train(NumEpoch);
    VectorND<T> acc_valid(NumEpoch);
    for (size_t epoch = 0; epoch < NumEpoch; ++epoch) {
        const auto nn_infer = MnistNet<T>(nn);
        loss_train[epoch] = nn_infer.loss(dataset.first);
        loss_valid[epoch] = nn_infer.loss(dataset.second);
        acc_train[epoch] = nn_infer.calcAccuracy(dataset.first);
        acc_valid[epoch] = nn_infer.calcAccuracy(dataset.second);

        nn.train_step_for<RandomSource, Sequential>(stepPerEpoch, BatchSize, dataset.first);
    }

    QApplication app(argc, argv);
    Plot* plot = new Plot(-0.5, NumEpoch - 0.5, 0, 3, 2, 1);
    auto& legend = plot->getLegend();
    legend.setAlignment(Qt::AlignRight);
    legend.setMarkerShape(QLegend::MarkerShapeFromSeries);
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
    axisRight->setMinorTickCount(4);
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
