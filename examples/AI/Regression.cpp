/*
 * Copyright 2022 WeiBo He.
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
#include "Physica/AI/Tensor.h"
#include "Physica/AI/RegressionDataset.h"
#include "Physica/AI/HyperParam/RandomSearch.h"
#include "Physica/AI/Model.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Utils/Random.h"

using namespace torch::data::datasets;
using namespace Physica::Core;
using namespace Physica::AI;

struct NetOptions;
class Net;

namespace Physica::AI::Internal {
    template<>
    class Traits<Net> {
    public:
        using DataSet = RegressionDataset;
        using HyperParams = NetOptions;
    };
}

struct NetOptions : public ParamSet<NetOptions> {
    size_t numEpoch;
    double lr;
    double gamma;

    friend std::ostream& operator<<(std::ostream& os, const NetOptions& options) {
        os << ' ' << options.numEpoch << ' ' << options.lr << ' ' << options.gamma << std::endl;
        return os;
    }

    template<class RandomGenerator>
    static NetOptions randomSet(RandomGenerator& gen) {
        NetOptions options{};
        {
            int arr[3]{300, 400, 500};
            std::uniform_int_distribution dist(0, 2);
            options.numEpoch = arr[dist(gen)];
        }
        {
            std::uniform_real_distribution<double> dist(-2, 1);
            options.lr = std::pow(10.0, dist(gen));
        }
        {
            std::uniform_real_distribution<double> dist(0.9, 1);
            options.gamma = dist(gen);
        }
        return options;
    }
};

class Net : public Model<Net> {
    using Base = Model<Net>;
    torch::nn::Linear fc1{nullptr}, fc2{nullptr}, fc3{nullptr};
public:
    Net(int inputs, int outputs) {
        fc1 = register_module("fc1", torch::nn::Linear(inputs, 256));
        fc2 = register_module("fc2", torch::nn::Linear(256, 256));
        fc3 = register_module("fc3", torch::nn::Linear(256, outputs));
    }

    void init() {
        fc1->weight.set_requires_grad(false).normal_(0, 0.01).set_requires_grad(true);
        fc2->weight.set_requires_grad(false).normal_(0, 0.01).set_requires_grad(true);
        fc3->weight.set_requires_grad(false).normal_(0, 0.01).set_requires_grad(true);
    }

    torch::Tensor forward(torch::Tensor x) {
        x = torch::relu(fc1->forward(x));
        x = torch::relu(fc2->forward(x));
        x = fc3->forward(x);
        return x;
    }

    void train(const RegressionDataset& dataset) {
        torch::data::DataLoaderOptions loader_option{};
        loader_option.batch_size() = 32;
        loader_option.enforce_ordering() = false;
        auto data_loader = torch::data::make_data_loader(
            const_cast<RegressionDataset&>(dataset).map(torch::data::transforms::Stack<>()), std::move(loader_option));

        torch::optim::AdamOptions adam_option{};
        adam_option.set_lr(active_params.lr);
        torch::optim::Adam optimizer(parameters(), std::move(adam_option));

        torch::optim::StepLR scheduler(optimizer, 1, active_params.gamma);

        for (size_t epoch = 1; epoch <= active_params.numEpoch; ++epoch) {
            for (auto& batch : *data_loader) {
                optimizer.zero_grad();
                torch::Tensor prediction = forward(batch.data);
                torch::Tensor mse = torch::mse_loss(prediction, batch.target);
                mse.backward();
                optimizer.step();
            }
            scheduler.step();
        }
    }

    double loss(const torch::Tensor& features, const torch::Tensor& labels) {
        const auto clipped_preds = forward(features).clamp(1);
        const auto rmse = torch::mse_loss(clipped_preds.log(), labels.log()).sqrt();
        return rmse.item().to<double>();
    }
};

RegressionDataset readTrainData() {
    using MatrixType = DenseMatrix<Scalar<Double, false>, DenseMatrixOption::Row | DenseMatrixOption::Vector>;
    MatrixType data(1460, 332);
    std::ifstream fin("../../data/train_num.csv");
    fin >> data;
    return RegressionDataset(toTensor(data.leftCols(331), at::kFloat),
                             toTensor(data.col(331).asVector(), at::kFloat).resize_({1460, 1}));
}

TensorDataset readTestData() {
    using MatrixType = DenseMatrix<Scalar<Double, false>, DenseMatrixOption::Row | DenseMatrixOption::Vector>;
    MatrixType data(1459, 331);
    std::ifstream fin("../../data/train_num.csv");
    fin >> data;
    return TensorDataset(toTensor(data));
}

int main() {
    auto net = std::make_shared<Net>(331, 1);
    KFold kFold(readTrainData(), 5);

    std::mt19937::result_type seed;
    Physica::Utils::Random::rdrand(seed);
    std::mt19937 gen(seed);

    NetOptions options{};
    options.numEpoch = 300;
    options.lr = 0.0212307;
    options.gamma = 0.95067;
    RandomSearch<Net> searcher(options, 0.158414);
    searcher.search(10, *net, kFold, gen);
    std::cout << searcher.getParams() << std::endl;
    std::cout << searcher.getScore() << std::endl;
    return 0;
}
