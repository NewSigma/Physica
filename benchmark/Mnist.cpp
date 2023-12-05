/*
 * Copyright 2023 WeiBo He. All rights reserved.
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
#include "Physica/Core/Math/Optimization/SGD.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"

using namespace Physica::Core;

template<class ScalarType> class MnistNet;

namespace Physica::Core::Internal {
    template<class T>
    class Traits<MnistNet<T>> {
    public:
        using ScalarType = T;
        using SampleType = Vector<T>;
        using LabelType = unsigned char;
    };
}

template<class ScalarType>
class MnistNet : public NetBase<MnistNet<ScalarType>> {
    using Base = NetBase<MnistNet<ScalarType>>;
    using typename Base::PlainScalar;
    using typename Base::VectorType;
    using typename Base::SampleType;
    using typename Base::LabelType;

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
    [[nodiscard]] VectorType forward(const VectorType& x) const {
        VectorType result = relu(layer1.forward(x));
        result = relu(layer2.forward(result));
        result = layer3.forward(result);
        return result;
    }

    template<class Optimizer>
    void opt_step(const Optimizer& opt) {
        layer1.opt_step(opt);
        layer2.opt_step(opt);
        layer3.opt_step(opt);
    }

    [[nodiscard]] ScalarType loss(const SampleType& sample, LabelType label) const {
        return Loss<ScalarType>::crossEntropy(forward(sample), label);
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
        const auto& testSamples = dataset.getTestDatas();
        const auto& testLabels = dataset.getTestLabels();
        const size_t numTestData = dataset.getNumTestData();
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
using Optimizer = SGD<PlainScalar>;
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

    auto& gen = RandomPoolType::getGen();
    const auto dataset = makeDataset();
    auto nn = MnistNet<ScalarType>(512, 512, gen);
    auto opt = Optimizer(0.01, batchSize);

    nn.train_step<Dataset, Optimizer, RandomPoolType, SequentialExecutor>(dataset, opt);
    return 0;
}
