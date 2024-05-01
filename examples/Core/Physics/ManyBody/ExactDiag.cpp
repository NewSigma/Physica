/*
 * Copyright 2024 WeiBo He.
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
#include <QApplication>
#include "Physica/Core/Math/Random/RandomPool.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/JacobiDavidson.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/SymmEigenSolver.h"
#include "Physica/Core/Physics/ManyBody/Hubbard1D.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;
using namespace Physica::Gui;
using ScalarType = Scalar<Double>;
using VectorType = Vector<ScalarType>;
using RandomPoolType = RandomPool<std::mt19937>;
constexpr unsigned int NumSite = 8;
constexpr unsigned int NumSpinUp = NumSite / 2;
constexpr unsigned int NumSpinDown = NumSite / 2;
constexpr double HoppingT = 1.0;

namespace Physica::Core {
    class HubbardS;

    namespace Internal {
        template<>
        class Traits<HubbardS> : public Traits<Hubbard1D<ScalarType>> {};
    }

    class HubbardS : public LatticeHamilton<HubbardS> {
        using Base = LatticeHamilton<HubbardS>;
        using typename Base::ScalarType;
        
        Hubbard1D<ScalarType> impl;
    public:
        HubbardS(Hubbard1D<ScalarType> impl_) : Base(impl_.getLattice()) {
            impl = std::move(impl_);
        }
        /* Operators */
        template<class VectorType>
        Vector<ScalarType> operator*(const RValueVector<VectorType>& v) const {
            const size_t length = v.getLength();
            const auto numSite = impl.getNumSuperCellSite();
            assert(Base::getColumn() == length && "[Error]: Dimensions do not match");
            Vector<ScalarType> result(length, 0);
            for (size_t i = 0; i < length; ++i) {
                const ScalarType hoppingElem = -v.calc(i) * impl.getHoppingT();
                const auto state = indexToState(i);
                int numRepel = 0;
                for (unsigned int site = 0; site < numSite; ++site) {
                    const auto site1 = (site + 1) % numSite;
                    Base::stateAdd(result, state.hopUp(site, site1), hoppingElem);
                    Base::stateAdd(result, state.hopUp(site1, site), hoppingElem);
                    Base::stateAdd(result, state.hopDown(site, site1), hoppingElem);
                    Base::stateAdd(result, state.hopDown(site1, site), hoppingElem);
                    if (site % 2U == 0U)
                        numRepel += state.isUpOccupy(site) && state.isDownOccupy(site);
                }
                result[i] += v.calc(i) * impl.getRepelU() * ScalarType(numRepel);
            }
            return result;
        }
        /* Operations */
        [[nodiscard]] size_t stateToIndex(StateType state) const noexcept { return impl.stateToIndex(state); }
        [[nodiscard]] StateType indexToState(size_t index) const noexcept { return impl.indexToState(index); }
        /* Getters */
        [[nodiscard]] size_t getNumState() const noexcept { return impl.getNumState(); }
    };
}

int main(int argc, char** argv) {
    ThreadPool::numThreadRequired = 4;
    const VectorType repelUs = VectorType::linspace(0, 30, 40);
    VectorType localMoments0(repelUs.getLength());
    VectorType localMoments1(repelUs.getLength());
    ThreadExecutor::parallel_for([&repelUs, &localMoments0, &localMoments1](unsigned int i) {
        HubbardS model(Hubbard1D<ScalarType>({{NumSite}, 1}, HoppingT, repelUs[i], NumSpinUp, NumSpinDown));
        const size_t numState = model.getNumState();

        JacobiDavidson<ScalarType> jd(numState, 4);
        jd.compute(model, VectorType::random_uniform(numState, RandomPoolType::getInstance().getGen()));
        jd.sort();

        const auto col = jd.getEigenvectors().col(0);
        ScalarType localMoment0 = 0, localMoment1 = 0;
        for (size_t i = 0; i < col.getLength(); ++i) {
            const auto state = model.indexToState(i);
            localMoment0 += square(col[i] * ScalarType(state.isUpOccupy(0) - state.isDownOccupy(0)));
            localMoment1 += square(col[i] * ScalarType(state.isUpOccupy(1) - state.isDownOccupy(1)));
        }
        localMoments0[i] = localMoment0;
        localMoments1[i] = localMoment1;
    }, repelUs.getLength()).wait();

    QApplication app(argc, argv);
    Plot* plot = new Plot(-0.1, 30.1, 0.37, 0.58, 10, 0.05);
    auto* axisX = plot->getAxisX();
    auto* axisY = plot->getAxisY();
    axisX->setTitleText("U/t");
    axisX->setLabelFormat("%d");
    axisY->setTitleText("&lt;m<sub>i</sub><sup>2</sup>&gt;");
    axisY->setLabelFormat("%.2f");
    {
        auto& scatter = plot->scatter(repelUs, localMoments0);
        scatter.setColor(Qt::black);
        scatter.setMarkerSize(10);
        scatter.setName("U<sub>i</sub> = U");
    }
    {
        auto& scatter = plot->scatter(repelUs, localMoments1);
        scatter.setColor(Qt::blue);
        scatter.setMarkerSize(10);
        scatter.setMarkerShape(QScatterSeries::MarkerShapeRectangle);
        scatter.setName("U<sub>i</sub> = 0");
    }
    plot->show();
    return QApplication::exec();
}
