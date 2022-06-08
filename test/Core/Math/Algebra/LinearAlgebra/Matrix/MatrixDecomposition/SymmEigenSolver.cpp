#include "Physica/Utils/TestHelper.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Transpose.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/SymmEigenSolver.h"

using namespace Physica::Core;

template<class MatrixType>
bool eigenTest(const MatrixType& mat, double precision) {
    using ScalarType = typename MatrixType::ScalarType;
    using RealType = typename ScalarType::RealType;
    using ComplexVector = Vector<ComplexScalar<RealType>, MatrixType::RowAtCompile, MatrixType::MaxRowAtCompile>;
    using EigenvectorMatrix = typename SymmEigenSolver<MatrixType>::EigenvectorMatrix;

    SymmEigenSolver<MatrixType> solver = SymmEigenSolver<MatrixType>(mat, true);
    solver.sort();

    const size_t order = mat.getRow();
    EigenvectorMatrix eigenvectors = solver.getEigenvectors();
    for (size_t i = 0; i < order; ++i) {
        if (i > 1 && solver.getEigenvalues()[i - 1] > solver.getEigenvalues()[i])
            return false;
        ComplexVector v1 = mat * eigenvectors.col(i);
        ComplexVector v2 = solver.getEigenvalues()[i] * eigenvectors.col(i).asVector();
        if (!vectorNear(v1, v2, precision))
            return false;
    }
    return true;
}

int main() {
    using RealType = Scalar<Double, false>;
    {
        using MatrixType = DenseMatrix<RealType, DenseMatrixOption::Column | DenseMatrixOption::Vector, 3, 3>;
        const MatrixType mat1{{-0.590316, -2.19514, -2.37463},
                             {-1.25006, -0.297493, 1.40349},
                             {0.517063, -0.956614, -0.920775}};
        MatrixType mat = mat1 + mat1.transpose();
        if (!eigenTest(mat, 1E-14))
            return 1;
    }
}
