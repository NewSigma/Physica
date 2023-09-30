#include <fstream>
#include "Physica/Core/IO/Poscar.h"
#include "Physica/Core/Physics/Phonon/FrozenPhonon.h"
#include "Physica/Core/Physics/SolidState/BrillouInterpolate.h"
#include "Physica/Utils/Unix/TempFile.h"

using namespace Physica::Core;
using namespace Physica::Utils;
using ScalarType = Scalar<Double>;
using ComplexType = ComplexScalar<ScalarType>;
using VectorType = Vector<ScalarType>;
using Vector3D = Vector<ScalarType, 3>;
using PhononType = FrozenPhonon<ScalarType, ScalarType>;
using MDCellType = typename PhononType::MDCellType;
using Index3D = typename GridBase::Index3D;

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

int main() {
    Poscar<ScalarType> poscar{};
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
    {
        const VectorType data{0.01944055685, 0.01944982364, 0.01947054153, 0.01949142063, 0.01950011668, 0.01949142063, 0.01947054153, 0.01944982364, 0.02056772432, 0.02055863517, 0.02005045933, 0.01987060425, 0.01981818079, 0.01981866608, 0.01987191262, 0.02005321502, 0.02123840976, 0.02151070523, 0.02123340334, 0.02083860036, 0.0206245673, 0.02056524248, 0.02062685159, 0.02084260051, 0.02167299944, 0.02188787602, 0.02188538597, 0.02166774908, 0.0214625393, 0.02137046087, 0.02137231221, 0.02146722801, 0.02181802013, 0.02184351139, 0.02186800762, 0.02183939408, 0.02181225543, 0.02183939408, 0.02186800762, 0.02184351139, 0.02167299944, 0.02146722801, 0.02137231221, 0.02137046087, 0.0214625393, 0.02166774908, 0.02188538597, 0.02188787602, 0.02123840976, 0.02084260051, 0.02062685159, 0.02056524248, 0.0206245673, 0.02083860036, 0.02123340334, 0.02151070523, 0.02056772432, 0.02005321502, 0.01987191262, 0.01981866608, 0.01981818079, 0.01987060425, 0.02005045933, 0.02055863517};
        const Index3D gridDim = {8, 8, 1};
        RSpaceGrid<ComplexType> dataGrid(gridDim);
        dataGrid.flatten() = data;

        BrillouInterpolate<ComplexType> zone({16, 16, 1}, unitCell.getLattice(), 0, 1E6);
        zone.interpolate(dataGrid);

        const Vector<ScalarType> x1 = Vector<ScalarType>::linspace(0, 1, 40);
        for (size_t i = 0; i < x1.getLength(); ++i) {
            const auto value = zone(Vector3D{x1[i], 0, 0});
            if (!scalarNear(value.getImag(), ScalarType(0), 1E-5))
                return 1;
        }
    }
    return 0;
}
