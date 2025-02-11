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
#include <gperftools/profiler.h>
#include <QApplication>
#include "Physica/Core/Math/Statistics/Correlation.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/Thermostat/DoubleThermo.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Physics/MD/ForceModel/LJModel1.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using KineticModel = FreeModel<ScalarType, 3, 1, RPMDIntegrator::Exact>;
using ThermoType = DoubleThermo<KineticModel>;
using MDType = RPMD<ScalarType, 3, 1>;
using RandomType = Random<MT19937>;
constexpr double thermostatTime = PhyConst<AU>::secondToTime(0.1 * 1E-12);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-14);
constexpr double latticeConst = PhyConst<AU>::angstormToBohr(5.67);
constexpr double temperatureT = PhyConst<AU>::kToTemperature(186);
constexpr size_t cellSize = 3;
constexpr size_t numReplica = 1;
constexpr size_t numSample = 300;
constexpr size_t numStep = 800;
constexpr size_t numSystem = 8;

namespace Physica {
    namespace Core {
        class ForceModel : public LJModel1<ScalarType> {
            constexpr static double lj_sigma = PhyConst<AU>::angstormToBohr(4.13);
            constexpr static double lj_epsilon = PhyConst<AU>::kToTemperature(165.9);
            constexpr static double pair_cutoff = 15;

            using Base = LJModel1<ScalarType>;
        public:
            ForceModel() : Base(lj_sigma, lj_epsilon, pair_cutoff) {}
        };
    }

    template<>
    class Traits<ForceModel> : public Traits<LJModel1<ScalarType>> {};
}

MDCell<ScalarType> makeSystem() {
    using CrystalCellType = CrystalCell<ScalarType>;
    constexpr size_t MoleculePerCell = 4;

    CrystalCellType::LatticeMatrix lattice = CrystalCellType::LatticeMatrix::unitMatrix(3);
    lattice *= ScalarType(latticeConst);

    CrystalCellType::PositionMatrix pos(MoleculePerCell, 3);
    {
        auto atomPos = pos.row(0);
        atomPos[0] = 0;
        atomPos[1] = 0;
        atomPos[2] = 0;
    }
    {
        auto atomPos = pos.row(1);
        atomPos[0] = latticeConst * 0.5;
        atomPos[1] = latticeConst * 0.5;
        atomPos[2] = 0;
    }
    {
        auto atomPos = pos.row(2);
        atomPos[0] = latticeConst * 0.5;
        atomPos[1] = 0;
        atomPos[2] = latticeConst * 0.5;
    }
    {
        auto atomPos = pos.row(3);
        atomPos[0] = 0;
        atomPos[1] = latticeConst * 0.5;
        atomPos[2] = latticeConst * 0.5;
    }

    CrystalCellType::AtomicArray atomicNumbers(MoleculePerCell, 36);

    CrystalCellType cell({std::move(lattice), std::move(pos), CrystalCellType::Type::Cartesian}, std::move(atomicNumbers));
    cell.toSuperCell(cellSize, cellSize, cellSize);
    return MDCell<ScalarType>(std::move(cell));
}

int main(int argc, char** argv) {
    VectorType corr_fft(numStep, 0), corr_fft_var(numStep, 0), corr_dir(numStep, 0), corr_dir_var(numStep, 0);
    {
        ThreadPool::numThreadRequired = 4;
        std::mutex mutex{};
        size_t sys = 0;
        ThreadExecutor::parallel_for([&](unsigned int) {
            auto rpmd = MDType(makeSystem(), numReplica, numReplica, temperatureT, timeStep);
            rpmd.initMomentum<KineticModel, RandomType>();

            ForceModel forceModel{};
            KineticModel kineticModel(temperatureT, numReplica);
            ThermoType thermo(temperatureT, thermostatTime);

            auto& ringPolymer = rpmd.getRingPolymer();
            Correlation<ScalarType> sampler(numStep);
            VectorType corr_dir_sample(numStep, 0);
            for (size_t sample = 0; sample < numSample; ++sample) {
                rpmd.nvt_step_for<ThermoType, RandomType, KineticModel, ForceModel, SeqExecutor>(
                        PhyConst<AU>::secondToTime(2 * 1E-12), thermo, kineticModel, forceModel);
                const ScalarType p0 = ringPolymer.makeCentroidMomentum()(0, 0);
                for (unsigned int step = 0; step < numStep; ++step) {
                    const ScalarType p1 = ringPolymer.makeCentroidMomentum()(0, 0);
                    rpmd.nvt_step<ThermoType, RandomType, KineticModel, ForceModel, SeqExecutor>(
                            thermo, kineticModel, forceModel);
                    sampler.sample(p1);
                    toNextMean(corr_dir_sample[step], sample, p0 * p1);
                }
            }
            mutex.lock();
            toNextVariance(corr_fft_var, corr_fft, sys, sampler.makeCorr(false));
            toNextVariance(corr_dir_var, corr_dir, sys, corr_dir_sample);
            sys += 1;
            mutex.unlock();
        }, numSystem, ThreadPool::numThreadRequired).wait();
    }

    const VectorType t = VectorType::linspace(0, numStep / 100.0, numStep);
    QApplication app(argc, argv);
    {
        Plot* plot = new Plot(0, double(*t.crbegin()) + 0.01, -0.5, 1, 2, 0.5);
        plot->getChart()->legend()->setAlignment(Qt::AlignTop);
        plot->getChart()->legend()->setMarkerShape(QLegend::MarkerShapeFromSeries);
        auto* axisX = plot->getAxisX();
        auto* axisY = plot->getAxisY();
        axisX->setTitleText("t/ps");
        axisX->setLabelFormat("%d");
        axisY->setTitleText("C<sub>pp</sub>/arb. unit");
        axisY->setLabelFormat("%.1f");

        QColor color;
        {
            auto& line = plot->line(VectorType{0, plot->getMaxX()}, VectorType{0, 0});
            color = line.color();
            auto pen = line.pen();
            pen.setColor(Qt::black);
            pen.setStyle(Qt::DashLine);
            line.setPen(pen);
        }
        {
            const ScalarType factor = reciprocal(corr_dir[0]);
            corr_dir *= factor;
            const VectorType corr_dir_devia = sqrt(corr_dir_var) * factor;

            auto& area = plot->area_center(t, corr_dir, corr_dir_devia);
            QColor color1 = color;
            color1.setAlpha(75);
            area.setColor(color1);
            auto pen = area.pen();
            pen.setColor(color1);
            area.setPen(pen);

            auto& line = plot->line(t, corr_dir);
            line.setColor(color);
            line.setName("Direct");
            plot->show();
        }
        color = Qt::red;
        {
            const ScalarType factor = reciprocal(corr_fft[0]);
            corr_fft *= factor;
            const VectorType corr_fft_devia = sqrt(corr_fft_var) * factor;

            auto& area = plot->area_center(t, corr_fft, corr_fft_devia);
            QColor color1 = color;
            color1.setAlpha(75);
            area.setColor(color1);
            auto pen = area.pen();
            pen.setColor(color1);
            area.setPen(pen);

            auto& line = plot->line(t, corr_fft);
            line.setColor(color);
            line.setName("FFT");
            plot->show();
        }
    }
    return QApplication::exec();
}
