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
#include <random>
#include <iostream>
#include <QtWidgets/QApplication>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Statistics/NumCharacter.h"
#include "Physica/Utils/Random.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;
using namespace Physica::Gui;
using ScalarType = Scalar<Float, false>;
/**
 * Reference:
 * [1] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013:304-308
 */
class Ising {
    DenseMatrix<ScalarType> lattice;
    ScalarType couplingJ;
    ScalarType boltzmannK;
    ScalarType temperature;
    ScalarType energy;
public:
    Ising(uint64_t N, ScalarType couplingJ_, ScalarType boltzmannK_, ScalarType temperature_)
            : lattice(N, N)
            , couplingJ(couplingJ_)
            , boltzmannK(boltzmannK_)
            , temperature(temperature_)
            , energy(0) {}
    /* Operations */
    template<class RandomGenerator>
    void init(RandomGenerator& generator) {
        std::uniform_real_distribution<float> dist{};
        for (uint64_t i = 0; i < lattice.getColumn(); ++i)
            for (uint64_t j = 0; j < lattice.getRow(); ++j)
                lattice(j, i) = (dist(generator) > 0.5) ? 1 : -1;
        energy = 0;
    }

    template<class RandomGenerator>
    void step(uint64_t stepNum, RandomGenerator& generator) {
        uint64_t iteration = stepNum * lattice.getRow() * lattice.getColumn();

        std::uniform_int_distribution<size_t> int_dist(0, lattice.getRow() - 1);
        std::uniform_real_distribution<float> dist{};
        const ScalarType beta = reciprocal(boltzmannK * temperature);
        for (uint64_t _ = 0; _ < iteration; ++_) {
            const size_t i = int_dist(generator);
            const size_t j = int_dist(generator);

            const ScalarType deltaE = -deltaDotSpin(i, j) * couplingJ;
            if (!deltaE.isPositive() || dist(generator) < exp(-deltaE * beta)) {
                lattice(i, j) = -lattice(i, j);
                energy += deltaE;
            }
        }
    }
    /* Getters */
    [[nodiscard]] ScalarType meanSpin() const {
        return lattice.sum() / square(ScalarType(lattice.getRow()));
    }

    [[nodiscard]] ScalarType getEnergy() const {
        return energy;
    }
private:
    ScalarType deltaDotSpin(size_t i, size_t j) const {
        const size_t order_1 = lattice.getRow() - 1;
        ScalarType spin = lattice(i > 0 ? (i - 1) : order_1, j);
        spin += lattice(i < order_1 ? (i + 1) : 0, j);
        spin += lattice(i, j > 0 ? (j - 1) : order_1);
        spin += lattice(i, j < order_1 ? (j + 1) : 0);
        spin *= lattice(i, j);
        return -spin * 2;
    }
};

int main(int argc, char** argv) {
    std::mt19937::result_type seed;
    Physica::Utils::Random::rdrand(seed);
    std::mt19937 gen(seed);

    const int count = 30;
    Vector<ScalarType> t(count);
    Vector<ScalarType> Cv(count);
    {
        Vector<ScalarType> buffer(500);
        for (int i = 0; i < count; ++i) {
            t[i] = 1 + i * 0.2;
            Ising ising(20, 1, 1, t[i]);
            ising.init(gen);
            ising.step(2000, gen);
            for (size_t i = 0; i < 500; ++i) {
                ising.step(10, gen);
                buffer[i] = ising.getEnergy();
            }
            auto e_mean = mean(buffer);
            buffer = square(buffer - e_mean);
            Cv[i] = mean(buffer) / square(t[i]);
        }
    }
    QApplication app(argc, argv);
    Plot* plot = new Plot();
    plot->spline(t, Cv);
    plot->chart()->createDefaultAxes();
    plot->show();
    return QApplication::exec();
}
