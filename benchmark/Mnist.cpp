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
#include "Physica/Core/AI/NeuralNetwork/Layer/LinearLayer.h"
#include "Physica/Core/AI/NeuralNetwork/SimpleNet.h"
#include "Physica/Core/AI/NeuralNetwork/Loss.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Math/Optimization/Stochastic/SGD.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"

using namespace Physica::Core;

template<Scalar T> class MnistNet;

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

    LinearLayer<T> layer1;
    LinearLayer<T> layer2;
    LinearLayer<T> layer3;
public:
    MnistNet() = default;
    template<class RandomType>
    MnistNet(size_t width1, size_t width2, RandomType& gen)
            : layer1(Mnist::NumPixelInImage, width1)
            , layer2(width1, width2)
            , layer3(width2, 10) {
        auto dist = std::normal_distribution<float>(0, 0.01);
        layer1.random_any(dist, gen);
        layer2.random_any(dist, gen);
        layer3.random_any(dist, gen);
    }
    template<Scalar U>
    MnistNet(const MnistNet<U>& net) : layer1(net.layer1), layer2(net.layer2), layer3(net.layer3) {}
    MnistNet(const MnistNet& other) = default;
    MnistNet(MnistNet&&) noexcept = default;
    ~MnistNet() = default;
    /* Operators */
    MnistNet& operator=(MnistNet& obj) noexcept { swap(obj); return *this; }
    /* Operations */
    [[nodiscard]] OutputType forward(const InputType& x) const {
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

    [[nodiscard]] MnistNet copy() const {
        MnistNet result{};
        result.layer1 = layer1.copy();
        result.layer2 = layer2.copy();
        result.layer3 = layer3.copy();
        return result;
    }

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
    template<Scalar> friend class MnistNet;
};

using ValueType = float32;
using ScalarType = Diff<ValueType, DiffMode::Reverse, 1>;
using Dataset = typename Mnist::DatasetType<VectorND<ValueType>>;
using Optimizer = SGD<ScalarType>;
using RandomType = Random<std::mt19937>;
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

        auto& gen = RandomType::getInstance().getGen();
        const auto dataset = makeDataset();
        auto opt = Optimizer(0.01, batchSize);
        opt.recordBegin();
        auto nn = MnistNet<ScalarType>(512, 512, gen);
        opt.recordEnd();

        for (auto _ : state)
            nn.train_step<Dataset, Optimizer, RandomType, SequentialExecutor>(dataset, opt);
    }
}

BENCHMARK(main)->Name("Mnist")->Unit(benchmark::kSecond);
