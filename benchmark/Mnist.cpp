/*
 * Copyright 2023-2024 Weibo He.
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
#include <benchmark/benchmark.h>
#include <gperftools/profiler.h>
#include "Physica/Core/IO/Mnist.h"
#include "Physica/Core/ML/NeuralNetwork/Layer/LinearLayer.h"
#include "Physica/Core/ML/NeuralNetwork/SeqNet.h"
#include "Physica/Core/ML/NeuralNetwork/Loss.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Math/Optimization/Stochastic/SGD.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"

using namespace Physica;

template<Scalar T, RandomGenerator R> class MnistNet;

namespace Physica {
    template<class T, RandomGenerator R>
    class Traits<MnistNet<T, R>> : public Traits<LinearLayer<T>> {};
}

template<Scalar T, RandomGenerator R>
class MnistNet : public SeqNet<MnistNet<T, R>> {
    using Base = SeqNet<MnistNet<T, R>>;
    using typename Base::ValueType;
    using typename Base::OutputType;

    LinearLayer<T> layer1;
    LinearLayer<T> layer2;
    LinearLayer<T> layer3;
public:
    MnistNet() = default;
    MnistNet(size_t width1, size_t width2)
            : layer1(Mnist::NumPixelInImage, width1)
            , layer2(width1, width2)
            , layer3(width2, 10) {
        auto dist = std::normal_distribution<float>(0, 0.01);
        layer1.template random_any<decltype(dist), R>(dist);
        layer2.template random_any<decltype(dist), R>(dist);
        layer3.template random_any<decltype(dist), R>(dist);
    }
    template<Scalar U>
    MnistNet(const MnistNet<U, R>& net) : layer1(net.layer1), layer2(net.layer2), layer3(net.layer3) {}
    MnistNet(const MnistNet& other) = default;
    MnistNet(MnistNet&&) noexcept = default;
    ~MnistNet() = default;
    /* Operators */
    MnistNet& operator=(MnistNet& obj) noexcept { swap(obj); return *this; }
    /* Operations */
    template<Vector V>
    [[nodiscard]] OutputType forward(const V& x) const {
        OutputType result = relu(layer1.forward(x));
        result = relu(layer2.forward(result));
        result = layer3.forward(result);
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
        for (const auto& sample : testSamples)
            count += testLabels[Base::classify(VectorType(sample))] == 1;
        return ValueType(count) / ValueType(numTestData);
    }

    void swap(MnistNet& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        layer1.swap(obj.layer1);
        layer2.swap(obj.layer2);
        layer3.swap(obj.layer3);
    }
private:
    template<Scalar, RandomGenerator> friend class MnistNet;
};

using ValueType = float32;
using ScalarType = Diff<ValueType, DiffMode::Reverse, 1>;
using Dataset = Mnist::DatasetType<VectorND<ValueType>>;
using Optimizer = SGD<ScalarType>;
using RandomType = Random<MT19937>;
constexpr size_t batchSize = 9000;

namespace {
    Dataset makeDataset() {
        const Mnist mnist("/home/sigma/Documents/data");
        auto dataset = mnist.makeTrainDataset<VectorND<ValueType>>();
        for (size_t i = 0; i < dataset.getSize(); ++i) {
            auto& sample = dataset.getSamples()[i];
            sample = sample * ValueType(1.0 / 128) - ValueType(1);
        }
        return dataset;
    }

    static void main(benchmark::State& state) {
        ThreadPool::numThreadRequired = 4;
        Dataset dataset;
        try {
            dataset = makeDataset();
        }
        catch (std::exception& e) {
            state.SkipWithError(e.what());
        }
        auto opt = Optimizer(0.01, batchSize);
        opt.recordBegin();
        auto nn = MnistNet<ScalarType, RandomType>(512, 512);
        opt.recordEnd();

        for (auto _ : state)
            nn.train_step<Dataset, Optimizer, RandomType, SequentialExecutor>(dataset, opt);
    }
}

BENCHMARK(main)->Name("Mnist")->Unit(benchmark::kSecond);
