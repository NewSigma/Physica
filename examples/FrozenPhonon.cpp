#include <fstream>
#include <QApplication>
#include "Physica/Core/IO/Poscar.h"
#include "Physica/Core/Physics/Phonon/FrozenPhonon.h"
#include "Physica/Core/Physics/MD/ForceModel/Q_TIP4P.h"
#include "Physica/Gui/Plot/PhononPlot.h"
#include "Physica/Utils/Unix/TempFile.h"

using namespace Physica::Core;
using namespace Physica::Gui;
using namespace Physica::Utils;
using ScalarType = Scalar<Double>;
using VectorType = Vector<ScalarType>;
using Vector3D = Vector<ScalarType, 3>;
using PhononType = FrozenPhonon<ScalarType, ScalarType>;
using InterpolateMethod = typename PhononType::InterpolateMethod;
using MatrixType = typename PhononType::MatrixType;
using MDCellType = typename PhononType::MDCellType;
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
                           "3.822044611    1.62606883   3.092008352\n"
                           "3.822721243  0.6591201425   1.908767104\n"
                           "3.817448139   2.351599932 0.05880066007\n"
                           "3.816771507   3.318548203   1.242042065\n"
                           "5.337373257   4.210567951 0.01213156618\n"
                           "2.294668913   4.210659504 0.01174806431\n"
                           "4.593580723    3.98064661   3.138677597\n"
                           "3.043938398   3.980554581   3.139060974\n"
                           "3.823375463  0.6846598387   2.877068996\n"
                           "3.816117287   3.293009043  0.2737401426\n"
                           "3.820759535  0.5593996644   0.148500219\n"
                           "3.818733454   3.418268919   3.002308846";

class ForceModel : private Q_TIP4P<ScalarType, ScalarType> {
    using Base = Q_TIP4P<ScalarType, ScalarType>;
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
        base.updateLattice(cell.getLattice());
        return base.potentialEnergy_unsort(cell);
    }
    using Base::getNumMolecule;
};

int main(int argc, char** argv) {
    Poscar poscar{};
    {
        auto tmp = TempFile("/tmp/tmpXXXXXX");
        std::ofstream os(tmp.getName());
        os << data;
        os.close();

        std::ifstream is(tmp.getName());
        is >> poscar;
        is.close();
        poscar.normalize();
    }
    MDCellType unitCell(poscar);
    unitCell.scale(PhyConst<AU>::angstormToBohr(1));

    PhononType ph(unitCell, {2, 2, 1});
    MDCellType superCell = unitCell.makeSuperCell<ExtendCellOption::CellMajor>(ph.getSuperSize());
    
    ForceModel model(superCell, pair_cutoff);
    const auto fcMatrixGrid = ph.makeForceConstants(model, displace, 1E-10, 100);

    QApplication app(argc, argv);
    auto* plot = new PhononPlot<ScalarType, ScalarType>();
    plot->chart()->legend()->setVisible(false);
    plot->plotPathLine<InterpolateMethod::FFT>(ph, fcMatrixGrid, {0, 0, 0}, {0.5, 0, 0}, 40, "M");
    plot->plotPathLine<InterpolateMethod::FFT>(ph, fcMatrixGrid, {0.5, 0, 0}, {0.5, 0.5, 0}, 40, "X");
    plot->plotPathLine<InterpolateMethod::FFT>(ph, fcMatrixGrid, {0.5, 0.5, 0}, {0, 0, 0}, 40, " Γ ");
    plot->setMaxY(120);
    plot->setDeltaY(50);
    plot->show();
    return QApplication::exec();
}
