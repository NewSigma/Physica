/*
 * Copyright 2023-2025 Weibo He.
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
#include <benchmark/benchmark.h>
#include <gperftools/profiler.h>
#include "Physica/Core/IO/Mnist.h"
#include "Physica/Core/ML/NeuralNetwork/Layer/LinearLayer.h"
#include "Physica/Core/ML/NeuralNetwork/SeqNet.h"
#include "Physica/Core/ML/Optimizer/SGD.h"
#include "Physica/Core/Math/Random/Random.h"

using namespace Physica;
constexpr size_t BatchSize = 64;
constexpr int HiddenW = 512;
constexpr double LearnRate = 0.01;

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

        [[nodiscard]] size_t classify(const VectorND<Tv>& input) const {
            static_assert(!Base::IsTrain, "[Error]: It is suggested using eval mode to reduce memory use");
            const auto output = Base::forward(input);
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
using dfloat = Diff<T, DiffMode::Reverse, 1>;
using Dataset = Mnist::DatasetType<VectorND<T>>;
using RandomSource = Random<>;

namespace {
    Dataset makeDataset() {
        const Mnist mnist("/home/sigma/Documents/data");
        auto dataset = mnist.makeTrainDataset<VectorND<T>>();
        for (size_t i = 0; i < dataset.getSize(); ++i) {
            auto& sample = dataset.getSamples()[i];
            sample = sample * T(1.0 / 128) - T(1);
        }
        return dataset;
    }

    void bench(benchmark::State& state) {
        Dataset dataset;
        try {
            dataset = makeDataset();
        }
        catch (std::exception& e) {
            state.SkipWithError(e.what());
            return;
        }

        auto nn = MnistNet<dfloat>(RandomSource::getInstance());
        for (auto _ : state)
            nn.train_step<RandomSource, Sequential>(BatchSize, dataset);
    }
}

BENCHMARK(bench)->Name("Mnist");
