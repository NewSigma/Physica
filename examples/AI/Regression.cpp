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
#include <QtWidgets/QApplication>
#include <QtCharts/QValueAxis>
#include "Physica/Gui/Plot/Plot.h"

using namespace torch::data::datasets;
using namespace Physica::Core;
using namespace Physica::Gui;
using namespace Physica::AI;
constexpr size_t numFeature = 101;

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
    unsigned int numEpoch;
    int64_t layer_dim1;
    int64_t layer_dim2;
    double lr;
    unsigned int step;
    double gamma;

    friend std::istream& operator>>(std::istream& is, NetOptions& options) {
        is >> options.numEpoch
           >> options.layer_dim1
           >> options.layer_dim2
           >> options.lr
           >> options.step
           >> options.gamma;
        return is;
    }

    friend std::ostream& operator<<(std::ostream& os, const NetOptions& options) {
        os << ' ' << options.numEpoch
           << ' ' << options.layer_dim1
           << ' ' << options.layer_dim2
           << ' ' << options.lr
           << ' ' << options.step
           << ' ' << options.gamma << std::endl;
        return os;
    }

    template<class RandomGenerator>
    static NetOptions randomSet(RandomGenerator& gen) {
        NetOptions options{};
        {
            int64_t arr[]{16, 32, 64, 128, 256};
            std::uniform_int_distribution dist(0, 4);
            options.layer_dim1 = arr[dist(gen)];
            options.layer_dim2 = arr[dist(gen)];
        }
        {
            std::uniform_int_distribution dist(300, 700);
            options.numEpoch = dist(gen);
        }
        {
            std::uniform_real_distribution<double> dist(-2, 1);
            options.lr = std::pow(10.0, dist(gen));
        }
        {
            std::uniform_int_distribution dist(1, 20);
            options.step = dist(gen);
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
    torch::nn::Linear fc1, fc2, fc3;
    int inputs;
    int outputs;
public:
    Net(int inputs_, int outputs_) : fc1(nullptr), fc2(nullptr), fc3(nullptr), inputs(inputs_), outputs(outputs_) {}

    void init() {
        if (fc1.is_empty()) {
            fc1 = register_module("fc1", torch::nn::Linear(inputs, active_params.layer_dim1));
            fc2 = register_module("fc2", torch::nn::Linear(active_params.layer_dim1, active_params.layer_dim2));
            fc3 = register_module("fc3", torch::nn::Linear(active_params.layer_dim2, outputs));
        }
        else {
            *fc1 = torch::nn::LinearImpl(inputs, active_params.layer_dim1);
            *fc2 = torch::nn::LinearImpl(active_params.layer_dim1, active_params.layer_dim2);
            *fc3 = torch::nn::LinearImpl(active_params.layer_dim2, outputs);
        }
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

    void trainOn(const RegressionDataset& dataset) {
        if (fc1.is_empty())
            init();
        torch::data::DataLoaderOptions loader_option{};
        loader_option.batch_size() = 32;
        loader_option.enforce_ordering() = false;
        auto data_loader = torch::data::make_data_loader(
            const_cast<RegressionDataset&>(dataset).map(torch::data::transforms::Stack<>()), std::move(loader_option));

        torch::optim::AdamOptions adam_option{};
        adam_option.set_lr(active_params.lr);
        torch::optim::Adam optimizer(parameters(), std::move(adam_option));

        torch::optim::StepLR scheduler(optimizer, active_params.step, active_params.gamma);

        for (unsigned int epoch = 1; epoch <= active_params.numEpoch; ++epoch) {
            for (auto& batch : *data_loader) {
                optimizer.zero_grad();
                torch::Tensor prediction = forward(batch.data);
                torch::Tensor temp = (prediction - batch.target) / batch.target;
                torch::Tensor loss = temp.norm();
                loss.backward();
                optimizer.step();
            }
            scheduler.step();
        }
    }

    double loss(const torch::Tensor& features, const torch::Tensor& labels) {
        if (is_training())
            toEvalMode();
        const auto clipped_preds = eval(features).clamp(1);
        const auto rmse = torch::mse_loss(clipped_preds.log(), labels.log()).sqrt();
        return rmse.item().to<double>();
    }
};

RegressionDataset readTrainData() {
    using MatrixType = DenseMatrix<Scalar<Double, false>, MatrixOption::Row | MatrixOption::Vector>;
    MatrixType data(1460, numFeature + 1);
    std::ifstream fin("../../data/train_num.csv");
    fin >> data;
    return RegressionDataset(toTensor(data.leftCols(data.getColumn() - 1), at::kFloat),
                             toTensor(data.col(data.getColumn() - 1).asVector(), at::kFloat).resize_({1460, 1}));
}

TensorDataset readTestData() {
    using MatrixType = DenseMatrix<Scalar<Double, false>, MatrixOption::Row | MatrixOption::Vector>;
    MatrixType data(1459, numFeature);
    std::ifstream fin("../../data/test_num.csv");
    fin >> data;
    return TensorDataset(toTensor(data));
}

int main(int argc, char** argv) {
    Net net(numFeature, 1);
    auto& options = net.active_params;
    options.gamma = 0.99;
    options.step = 500;

    options.numEpoch = 25;
    options.lr = 0.002;
    options.layer_dim1 = 32;
    options.layer_dim2 = 64;
    KFold kFold(readTrainData(), 5);

    std::mt19937::result_type seed;
    Physica::Utils::Random::rdrand(seed);
    torch::manual_seed(seed);

    auto curve = kFold.makeLearningCurve(net, 500);

    QApplication app(argc, argv);
    QFont font;
    Plot* plot = new Plot();
    auto& chart = *plot->chart();
    chart.legend()->setVisible(false);
    {
        constexpr double minX = 0;
        constexpr double maxX = 10;
        constexpr double minY = 50;
        constexpr double maxY = 120;
        QValueAxis* axisX = new QValueAxis();
        font = axisX->labelsFont();
        font.setPointSize(15);
        axisX->setTickAnchor(0);
        axisX->setTickInterval(2);
        //axisX->setTickType(QValueAxis::TicksDynamic);
        axisX->setMinorGridLineVisible(false);
        axisX->setLinePenColor(Qt::black);
        axisX->setGridLineVisible(false);
        axisX->setLabelsFont(font);
        //axisX->setRange(minX, maxX);
        axisX->setTitleFont(font);
        QValueAxis* axisY = new QValueAxis();
        axisY->setTickAnchor(0);
        axisY->setTickInterval(10);
        //axisY->setTickType(QValueAxis::TicksDynamic);
        axisY->setMinorGridLineVisible(false);
        axisY->setMinorTickCount(4);
        axisY->setLinePenColor(Qt::black);
        axisY->setGridLineVisible(false);
        axisY->setMinorGridLineVisible(false);
        axisY->setLabelsFont(font);
        //axisY->setRange(minY, maxY);
        axisY->setTitleFont(font);

        chart.addAxis(axisX, Qt::AlignBottom);
        chart.addAxis(axisY, Qt::AlignLeft);

        {
            auto& scatter = plot->line(curve.first);
            scatter.attachAxis(axisX);
            scatter.attachAxis(axisY);
        }
        {
            auto& scatter = plot->line(curve.second);
            scatter.attachAxis(axisX);
            scatter.attachAxis(axisY);
        }
    }
    plot->show();
    return QApplication::exec();
}
