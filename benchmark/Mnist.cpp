/*
 * Copyright 2023-2024 WeiBo He. All rights reserved.
 *
 * This file is part of PhysicaNotes.
 */
#include <iostream>
#include <gperftools/profiler.h>
#include "Physica/Core/IO/Mnist.h"
#include "Physica/Core/AI/NeuralNetwork/Layer/LinearLayer.h"
#include "Physica/Core/AI/NeuralNetwork/NetBase.h"
#include "Physica/Core/AI/NeuralNetwork/Loss.h"
#include "Physica/Core/Math/Random/RandomPool.h"
#include "Physica/Core/Math/Optimization/Stochastic/SGD.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"

using namespace Physica::Core;

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

    LinearLayer<ScalarType> layer1;
    LinearLayer<ScalarType> layer2;
    LinearLayer<ScalarType> layer3;
public:
    MnistNet() = default;
    template<class RandomGenerator>
    MnistNet(size_t width1, size_t width2, RandomGenerator& gen)
            : layer1(Mnist::NumPixelInImage, width1)
            , layer2(width1, width2)
            , layer3(width2, 10) {
        auto dist = std::normal_distribution<float>(0, 0.01);
        layer1.random_any(dist, gen);
        layer2.random_any(dist, gen);
        layer3.random_any(dist, gen);
    }
    template<class OtherScalar>
    MnistNet(const MnistNet<OtherScalar>& net) : layer1(net.layer1), layer2(net.layer2), layer3(net.layer3) {}
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
    [[nodiscard]] ScalarType loss(const Dataset& dataset, size_t index) const {
        return Loss<ScalarType>::crossEntropy(forward(dataset.getSamples()[index]), dataset.getLabels()[index]);
    }

    template<class Dataset>
    [[nodiscard]] ScalarType loss(const Dataset& dataset) const { return Base::loss(dataset); }

    [[nodiscard]] MnistNet copy() const {
        MnistNet result{};
        result.layer1 = layer1.copy();
        result.layer2 = layer2.copy();
        result.layer3 = layer3.copy();
        return result;
    }

    template<class Dataset>
    PlainScalar calcAccuracy(const Dataset& dataset) const {
        const auto& testSamples = dataset.getSamples();
        const auto& testLabels = dataset.getLabels();
        const size_t numTestData = dataset.getSize();
        size_t count = 0;
        for (const auto& sample : testSamples)
            count += testLabels[Base::classify(VectorType(sample))] == 1;
        return PlainScalar(count) / PlainScalar(numTestData);
    }

    void swap(MnistNet& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Base::swap(obj);
        layer1.swap(obj.layer1);
        layer2.swap(obj.layer2);
        layer3.swap(obj.layer3);
    }
private:
    template<class OtherScalar> friend class MnistNet;
};

using PlainScalar = Scalar<Float>;
using ScalarType = Differentiable<PlainScalar, DiffMode::Reverse>;
using Dataset = typename Mnist::DatasetType<PlainScalar>;
using Optimizer = SGD<ScalarType>;
using RandomPoolType = RandomPool<std::mt19937>;
constexpr size_t numEpoch = 10;
constexpr size_t batchSize = 9000;

Dataset makeDataset() {
    const Mnist mnist("/home/sigma/Documents/data");
    auto dataset = mnist.makeTrainDataset<PlainScalar>();
    for (size_t i = 0; i < dataset.getSize(); ++i) {
        auto& sample = dataset.getSamples()[i];
        sample = sample * PlainScalar(1.0 / 128) - PlainScalar(1);
    }
    return dataset;
}

int main() {
    ThreadPool::numThreadRequired = 4;

    auto& gen = RandomPoolType::getInstance().getGen();
    const auto dataset = makeDataset();
    auto opt = Optimizer(0.01, batchSize);
    opt.recordBegin();
    auto nn = MnistNet<ScalarType>(512, 512, gen);
    opt.recordEnd();

    nn.train_step<Dataset, Optimizer, RandomPoolType, SequentialExecutor>(dataset, opt);
    return 0;
}
