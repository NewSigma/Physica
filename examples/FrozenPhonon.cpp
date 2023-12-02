#include <fstream>
#include <QApplication>
#include "Physica/Core/IO/Poscar.h"
#include "Physica/Core/Physics/MD/ForceModel/Q_TIP4P.h"
#include "Physica/Core/Physics/MD/ForceModel/Ewald/Ewald.h"
#include "Physica/Core/Physics/Phonon/FrozenPhonon.h"
#include "Physica/Core/Physics/Phonon/PhononDOS.h"
#include "Physica/Gui/Plot/PhononPlot.h"
#include "Physica/Utils/Unix/TempFile.h"

using namespace Physica::Core;
using namespace Physica::Gui;
using namespace Physica::Utils;
using ScalarType = Scalar<Double>;
using ComplexType = ComplexScalar<ScalarType>;
using VectorType = Vector<ScalarType>;
using Vector3D = Vector<ScalarType, 3>;
using PhononType = FrozenPhonon<ScalarType>;
using MatrixType = typename PhononType::MatrixType;
using MDCellType = typename PhononType::MDCellType;
using Index3D = typename GridBase::Index3D;
constexpr double pair_cutoff = PhyConst<AU>::angstormToBohr(9);
constexpr double displace = PhyConst<AU>::angstormToBohr(0.001);

const static char* data = "Structure\n"
                           "1\n"
                           "4.591253757           0           0\n"
                           "2.290947199  4.21609211           0\n"
                           "1.993453026 1.389889121          18\n"
                           "H O\n"
                           "8 4\n"
                           "Cartesian\n"
                           "3.82392771      1.627786996     3.097418093\n"
                           "3.825170743     0.6609370111    1.91455944\n"
                           "3.821134627     2.355310853     0.05898047382\n"
                           "3.819906378     3.322169175     1.24183939\n"
                           "5.340537394     4.216015464     0.01506002955\n"
                           "2.298689371     4.216120568     0.01476307969\n"
                           "4.595492317     3.983179478     3.141328029\n"
                           "3.046108593     3.983091746     3.141638832\n"
                           "3.825339204     0.6863512929    2.882796309\n"
                           "3.819702744     3.296745616     0.2735997228\n"
                           "3.824291156     0.562119077     0.152752783\n"
                           "3.820790917     3.420982462     3.003610409";

class ForceModel : private Q_TIP4P<ScalarType, Ewald<ScalarType>> {
    using Base = Q_TIP4P<ScalarType, Ewald<ScalarType>>;
    using typename Base::MDCellType;
public:
    using Base::Base;
    /* Operations */
    template<class Executor, bool IsSmallCell = true>
    [[nodiscard]] Vector<ScalarType> force(const MDCellType& cell) const {
        static_assert(IsSmallCell);
        return Base::force_unsort<Executor, true>(cell);
    }
    ScalarType potentialEnergy(const MDCellType& cell) const {
        Base base(*this);
        base.setLattice(cell.getLattice());
        return base.potentialEnergy_unsort(cell);
    }
    using Base::getNumMolecule;
};

int main(int argc, char** argv) {
    Poscar<ScalarType> poscar{};
    try {
        auto tmp = TempFile("/tmp/tmpXXXXXX");
        std::ofstream os(tmp.getName());
        os << data;
        os.close();

        std::ifstream is(tmp.getName());
        is >> poscar;
        is.close();
        poscar.normalize();
    }
    catch (std::exception& e) {
        std::cout << e.what() << std::endl;
        return 1;
    }
    MDCellType unitCell(poscar);
    unitCell.scale(PhyConst<AU>::angstormToBohr(1));

    PhononType ph(unitCell, {8, 8, 1}, displace);
    MDCellType superCell = unitCell.makeSuperCell<ExtendCellOption::CellMajor>(ph.getSuperSize());
    
    ForceModel model(superCell, pair_cutoff, {});
    auto fcMatrixGrid = ph.makeForceConstants(model);
    ph.applyTranslate(fcMatrixGrid, 1E-10, 100);

    QApplication app(argc, argv);
    PhononPlot<ScalarType>* phPlot;
    Plot* dosPlot;
    {
        phPlot = new PhononPlot<ScalarType>();
        phPlot->chart()->legend()->setVisible(false);
        phPlot->plotPathLine(ph, fcMatrixGrid, {0, 0, 0}, {0.5, 0, 0}, 40, "M");
        phPlot->plotPathLine(ph, fcMatrixGrid, {0.5, 0, 0}, {0.5, 0.5, 0}, 40, "X");
        phPlot->plotPathLine(ph, fcMatrixGrid, {0.5, 0.5, 0}, {0, 0, 0}, 40, " Γ ");
        phPlot->setMaxY(120);
        phPlot->setDeltaY(50);
    }
    {
        PhononDOS<ScalarType> dosSolver(ph.getUnitCell(), ph.getSuperSize(), fcMatrixGrid, {32, 32, 1});
        auto freq = Vector<ScalarType>::linspace(0, 110, 300);
        freq *= ScalarType(PhyConst<AU>::THzToFreq(1));
        Vector<ScalarType> dos(freq.getLength());
        for (size_t i = 0; i < dos.getLength(); ++i)
            dos[i] = dosSolver.calcDOS(freq[i]);
        dos.toUnit();

        dosPlot = new Plot(0, 0.4, 0, 120, 0.1, 50);
        dosPlot->chart()->legend()->setVisible(false);
        auto* axisX = dosPlot->getAxisX();
        auto* axisY = dosPlot->getAxisY();
        axisX->setTitleText("Normalized DOS/THz<sup>-1</sup>");
        axisX->setLabelFormat("%.1f");
        axisY->setTitleText("Frequency/THz");
        axisY->setLabelFormat("%d");

        freq *= ScalarType(PhyConst<AU>::freqToTHz(1));
        dosPlot->line(dos, freq).setColor(Qt::black);
    }
    phPlot->show();
    dosPlot->show();
    return QApplication::exec();
}
