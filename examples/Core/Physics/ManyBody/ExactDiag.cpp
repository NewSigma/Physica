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
#include <Physica/Core/Math/Random/RandomPool.h>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Eigen/JacobiDavidson.h>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Eigen/SymmEigenSolver.h>
#include <Physica/Core/Physics/ManyBody/Hubbard.h>
#include <Physica/Core/Parallel/Executor/ThreadExecutor.h>
#include <Physica/Gui/Plot/Plot.h>

using namespace Physica::Core;
using namespace Physica::Gui;
using ScalarType = Scalar<Double>;
using VectorType = Vector<ScalarType>;
using RandomPoolType = RandomPool<std::mt19937>;
constexpr unsigned int NumSite = 8;
constexpr unsigned int NumSpinUp = NumSite / 2;
constexpr unsigned int NumSpinDown = NumSite / 2;
constexpr double HoppingT = 1.0;
using ReprType = SpinRepr<NumSpinUp == NumSpinDown>;
using HubbardType = Hubbard<ScalarType, ReprType>;

namespace Physica::Core {
    class Hamilton;

    namespace Internal {
        template<>
        class Traits<Hamilton> : public Traits<Hubbard<ScalarType, ReprType>> {};
    }

    class Hamilton : public RValueMatrix<Hamilton>, private HubbardType {
        using Base = RValueMatrix<Hamilton>;
        using ModelBase = HubbardType;
        using typename Base::ScalarType;

        ScalarType magneticH;
    public:
        Hamilton(HubbardType base, ScalarType magneticH_) : ModelBase(std::move(base)), magneticH(magneticH_) {}
        /* Operators */
        template<class VectorType>
        Vector<ScalarType> operator*(const RValueVector<VectorType>& v) const {
            const size_t length = v.getLength();
            const auto numSite = ModelBase::getNumSuperCellSite();
            assert(getColumn() == length && "[Error]: Dimensions do not match");
            Vector<ScalarType> result(length, 0);
            for (size_t i = 0; i < length; ++i) {
                const ScalarType hoppingElem = -v.calc(i) * ModelBase::getHoppingT();
                const auto state = getRepr()[i];
                int numRepel = 0;
                int numPolar = 0;
                for (unsigned int site = 0; site < numSite; ++site) {
                    const auto site1 = (site + 1) % numSite;
                    sumHopping(result, hoppingElem, state.hopUp(site, site1));
                    sumHopping(result, hoppingElem, state.hopUp(site1, site));
                    sumHopping(result, hoppingElem, state.hopDown(site, site1));
                    sumHopping(result, hoppingElem, state.hopDown(site1, site));
                    const auto upOccupy = state.isUpOccupy(site);
                    const auto downOccupy = state.isDownOccupy(site);
                    numRepel += (site % 2U == 0U) && upOccupy && downOccupy;
                    numPolar += int(upOccupy) - int(downOccupy);
                }
                result[i] += v.calc(i) * (ModelBase::getRepelU() * ScalarType(numRepel) - magneticH * ScalarType(numPolar));
            }
            return result;
        }
        /* Operations */
        using Base::transpose;
        [[nodiscard]] const Hamilton& hermite() const noexcept { return *this; }
        /* Getters */
        using ModelBase::getRow;
        using ModelBase::getColumn;
        using ModelBase::getRepr;
        using ModelBase::getNumState;
    };
}

int main(int argc, char** argv) {
    ThreadPool::numThreadRequired = 4;
    const VectorType repelUs = VectorType::linspace(0, 30, 40);
    VectorType localMoments0(repelUs.getLength());
    VectorType localMoments1(repelUs.getLength());
    ThreadExecutor::parallel_for([&repelUs, &localMoments0, &localMoments1](unsigned int i) {
        ReprType repr(NumSite, NumSpinUp, NumSpinDown);
        Hamilton model(HubbardType({{NumSite}, 1}, std::move(repr), HoppingT, repelUs[i]), ScalarType(0));

        const size_t numState = model.getNumState();
        JacobiDavidson<ScalarType> jd(numState, 4);
        jd.template compute<Hamilton>(model, VectorType::random_uniform(numState, RandomPoolType::getInstance().getGen()));
        jd.sort();

        const auto col = jd.getEigenvectors().col(0);
        ScalarType localMoment0 = 0, localMoment1 = 0;
        for (size_t i = 0; i < col.getLength(); ++i) {
            const auto state = model.getRepr()[i];
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
