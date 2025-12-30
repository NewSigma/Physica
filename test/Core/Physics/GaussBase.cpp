/*
 * Copyright 2021-2025 Weibo He.
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
#include "Physica/Core/Physics/ElectronicStructure/HF/GaussBase.h"
#include "Test.h"

using namespace Physica;

using ScalarType = float64;

namespace Physica {
    class Test {
    public:
        [[nodiscard]] static ScalarType helper_F(size_t v, const ScalarType& t) { return GaussBase<ScalarType>::helper_F(v, t); }
    };
}

namespace {
    void test_helper_F() {
        expect(scalarNear(Test::helper_F(0, 39.03714925), ScalarType(0.1418423419), 1E-10));
    }

    ScalarType overlap_1s_1s(const ScalarType& alpha1,
                            const Vector3D<ScalarType>& v1,
                            const ScalarType& alpha2,
                            const Vector3D<ScalarType>& v2) {
        const ScalarType alpha_sum = alpha1 + alpha2;
        const ScalarType temp = ScalarType(M_PI) / alpha_sum;
        const ScalarType factor = temp * sqrt(temp);
        return factor * exp(-alpha1 * alpha2 / alpha_sum * (v1 - v2).squaredNorm());
    }

    ScalarType kinetic_1s_1s(const ScalarType& alpha1,
                            const Vector3D<ScalarType>& v1,
                            const ScalarType& alpha2,
                            const Vector3D<ScalarType>& v2) {
        const ScalarType alpha_sum = alpha1 + alpha2;
        const ScalarType temp = ScalarType(M_PI) / alpha_sum;
        const ScalarType factor = temp * sqrt(temp);
        const ScalarType temp1 = alpha1 * alpha2 / alpha_sum;
        const ScalarType squaredNorm = (v1 - v2).squaredNorm();
        return factor * exp(-temp1 * squaredNorm) * temp1 * (ScalarType(6) - ScalarType(4) * temp1 * squaredNorm) * ScalarType(0.5);
    }

    ScalarType attraction_1s_1s(const ScalarType& alpha1,
                            const Vector3D<ScalarType>& v1,
                            const ScalarType& alpha2,
                            const Vector3D<ScalarType>& v2,
                            const Vector3D<ScalarType>& corePos) {
        const ScalarType alpha_sum = alpha1 + alpha2;
        const ScalarType factor = ScalarType(2 * M_PI) / alpha_sum;
        const ScalarType temp1 = alpha1 / alpha_sum;
        const ScalarType temp2 = ScalarType(1) - temp1;
        const Vector3D<ScalarType> vector_p = temp1 * v1 + temp2 * v2;
        const ScalarType squaredNorm = (v1 - v2).squaredNorm();
        return factor * exp(-temp1 * alpha2 * squaredNorm) * Test::helper_F(0, alpha_sum * (vector_p - corePos).squaredNorm());
    }

    ScalarType repulsion_1s_1s(const ScalarType& alpha1,
                            const Vector3D<ScalarType>& v1,
                            const ScalarType& alpha2,
                            const Vector3D<ScalarType>& v2,
                            const ScalarType& alpha3,
                            const Vector3D<ScalarType>& v3,
                            const ScalarType& alpha4,
                            const Vector3D<ScalarType>& v4) {
        const ScalarType alpha_sum1 = alpha1 + alpha3;
        const ScalarType alpha_sum2 = alpha2 + alpha4;
        const ScalarType factor = ScalarType(2 * std::pow(M_PI, 2.5)) / (alpha_sum1 * alpha_sum2 * sqrt(alpha_sum1 + alpha_sum2));
        const ScalarType temp1 = alpha1 / alpha_sum1;
        const ScalarType temp2 = ScalarType(1) - temp1;
        const Vector3D<ScalarType> vector_p = temp1 * v1 + temp2 * v3;
        const ScalarType squaredNorm1 = (v1 - v3).squaredNorm();
        const ScalarType temp3 = alpha2 / alpha_sum2;
        const ScalarType temp4 = ScalarType(1) - temp3;
        const Vector3D<ScalarType> vector_q = temp3 * v2 + temp4 * v4;
        const ScalarType squaredNorm2 = (v2 - v4).squaredNorm();
        return factor * exp(-temp1 * alpha3 * squaredNorm1) * exp(-temp3 * alpha4 * squaredNorm2)
            * Test::helper_F(0, alpha_sum1 * alpha_sum2 * (vector_p - vector_q).squaredNorm() / (alpha_sum1 + alpha_sum2));
    }
}

int main() {
    using BaseFunc = GaussBase<ScalarType>;
    test_helper_F();

    {
        auto alpha1 = ScalarType(1.25);
        Vector3D<ScalarType> v1{2, 5, -1};
        GaussBase<ScalarType> base1 = GaussBase<ScalarType>(v1, alpha1, 0, 0, 0);
        auto alpha2 = ScalarType(0.76);
        Vector3D<ScalarType> v2{-3, 6, 1};
        GaussBase<ScalarType> base2 = GaussBase<ScalarType>(v2, alpha2, 0, 0, 0);
        expect(scalarNear(BaseFunc::overlap(base1, base2), overlap_1s_1s(alpha1, v1, alpha2, v2), 1E-14));
        expect(scalarNear(BaseFunc::kinetic(base1, base2), kinetic_1s_1s(alpha1, v1, alpha2, v2), 1E-14));

        Vector3D<ScalarType> v3{1.5, 1.7, -0.4};
        expect(scalarNear(BaseFunc::nuclearAttraction(base1, base2, v3), attraction_1s_1s(alpha1, v1, alpha2, v2, v3), 1E-14));

        auto alpha3 = ScalarType(3.78);
        auto alpha4 = ScalarType(11.7);
        Vector3D<ScalarType> v4{2.7, 0, -3};
        GaussBase<ScalarType> base3 = GaussBase<ScalarType>(v3, alpha3, 0, 0, 0);
        GaussBase<ScalarType> base4 = GaussBase<ScalarType>(v4, alpha4, 0, 0, 0);
        expect(scalarNear(BaseFunc::electronRepulsion(base1, base2, base3, base4),
                          repulsion_1s_1s(alpha1, v1, alpha2, v2, alpha3, v3, alpha4, v4),
                          1E-14));
    }
    {
        Vector3D<ScalarType> v{0, 0, 0};
        GaussBase<ScalarType> base1(v, 0.89, 1, 0, 0);
        GaussBase<ScalarType> base2(v, 12.7, 2, 0, 0);
        expect(scalarNear(BaseFunc::overlap(base1, base2), ScalarType(0), 1E-15));

        GaussBase<ScalarType> base3(v, 0.89, 2, 0, 0);
        expect(scalarNear(BaseFunc::overlap(base2, base3), ScalarType(0.0004513547048841694), 1E-15));
    }
    {
        Vector3D<ScalarType> v{0, 0, 0};
        GaussBase<ScalarType> base1(v, 0.298073, 0, 0, 0);
        GaussBase<ScalarType> base2(v, 0.298073, 0, 0, 0);
        expect(scalarNear(BaseFunc::nuclearAttraction(base1, base2, v), ScalarType(10.539675360028560918), 1E-15));
    }
    {
        const auto alpha1 = ScalarType(13.00773);
        const auto alpha2 = ScalarType(0.444529);
        const Vector3D<ScalarType> v1{0, 0, 0};
        const Vector3D<ScalarType> v2{0, 0, 1};
        const GaussBase<ScalarType> base1(v1, alpha1, 0, 0, 0);
        const GaussBase<ScalarType> base2(v2, alpha2, 0, 0, 0);
        expect(scalarNear(BaseFunc::electronRepulsion(base1, base2, base2, base1),
                          repulsion_1s_1s(alpha1, v1, alpha2, v2, alpha2, v2, alpha1, v1),
                          1E-14));
    }
    return 0;
}
