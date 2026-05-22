/*
 * Copyright 2022-2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica;
using T = float32;
using RandomSource = Random<>;
constexpr int NumPoint = 120;
constexpr int NumSample = 5000;

namespace {
    class Ising final {
        DenseMatrix<T> lattice;
        T energy;
        T betaJ;
    public:
        Ising(uint64_t N, T betaJ_) : lattice(N, N), energy(0), betaJ(betaJ_) {}
        /* Operations */
        template<RNG R>
        void init() {
            std::uniform_real_distribution<float> dist{};
            for (auto& spin : lattice.flatten())
                spin = (dist(R::getInstance()) > 0.5) ? 1 : -1;
            energy = 0;
        }

        template<RNG R>
        void step(uint64_t stepNum) {
            uint64_t iteration = stepNum * lattice.getRow() * lattice.getCol();

            std::uniform_int_distribution<size_t> int_dist(0, lattice.getRow() - 1);
            auto& gen = R::getInstance();
            for (uint64_t _ = 0; _ < iteration; ++_) {
                const size_t i = int_dist(gen);
                const size_t j = int_dist(gen);

                const T deltaE = -deltaDotSpin(i, j);
                if (!deltaE.isPositive() || T::random_uniform<R>() < exp(-deltaE * betaJ)) {
                    lattice[i, j] = -lattice[i, j];
                    energy += deltaE;
                }
            }
        }
        /* Getters */
        [[nodiscard]] T meanSpin() const {
            return lattice.sum() / square(T(lattice.getRow()));
        }

        [[nodiscard]] T getEnergy() const {
            return energy;
        }

        [[nodiscard]] size_t getNumSite() const noexcept {
            return lattice.getSize();
        }
    private:
        [[nodiscard]] T deltaDotSpin(size_t i, size_t j) const {
            const size_t order_1 = lattice.getRow() - 1;
            T spin = lattice[i > 0 ? (i - 1) : order_1, j];
            spin += lattice[i < order_1 ? (i + 1) : 0, j];
            spin += lattice[i, j > 0 ? (j - 1) : order_1];
            spin += lattice[i, j < order_1 ? (j + 1) : 0];
            spin *= lattice[i, j];
            return -spin * 2;
        }
    };
}
/**
 * Reference:
 * [1] J. H. Thijssen. Computational Physics[M]. London: Cambridge University Press, 2013:304-308
 */
int main(int argc, char** argv) {
    ThreadPool::numThreadRequired = 4;
    const auto t = VectorND<T>::linspace(1, 7, NumPoint);
    VectorND<T> Cv(NumPoint);
    parallel_for<Thread>([&](size_t i) {
        Ising ising(20, reciprocal(t[i]));
        ising.init<RandomSource>();
        ising.step<RandomSource>(2000);

        T energy = 0;
        T energy2 = 0;
        for (size_t i = 0; i < NumSample; ++i) {
            ising.step<RandomSource>(10);
            energy.toNextMean(i, ising.getEnergy());
            energy2.toNextMean(i, square(ising.getEnergy()));
        }
        Cv[i] = (energy2 - square(energy)) / (square(t[i]) * T(ising.getNumSite()));
    }, NumPoint, ThreadPool::numThreadRequired).wait();

    QApplication app(argc, argv);
    Plot* plot = new Plot(1, 7, 0, 2, 1, 1);
    plot->getLegend().hide();
    plot->getAxisX()->setLabelFormat("%d");
    plot->getAxisY()->setLabelFormat("%d");
    plot->getAxisX()->setTitleText("T");
    plot->getAxisY()->setTitleText("C<sub>v</sub>");
    plot->line(t, Cv);
    plot->show();
    return QApplication::exec();
}
