/*
 * Copyright 2026 Weibo He.
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
#pragma once

#include <iostream>
#include "Vector/DenseVector.h"
#include "LinearSystem/GMRES.h"
#include "Physica/Core/Math/Calculus/Integrate/Integrate.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiagMatrix.h"

namespace Physica {
    /**
     * \class SparseDet: lnAbsDet for general large scale sparse matrix
     *
     * Reformulation of the lnAbsDet by an integral representation and stochastic expressions, refer to [1]
     *
     * Reference:
     * [1] Phys. Rev. E 92, 013302 (2015); https://doi.org/10.1103/PhysRevE.92.013302
     */
    template<Scalar T>
    class SparseDet {
        using This = SparseDet<T>;
        using Tr = T::RealType;
        using Trv = Tr::ValueType;

        GMRES<T> gmres;
        DiagMatrix<T> precond;
        DiagMatrix<T> matD;
        VectorND<T> rand;
        VectorND<T> solX;
    public:
        SparseDet(GMRES<T> gmres_);
        SparseDet(const This&) = default;
        SparseDet(This&&) noexcept = default;
        ~SparseDet() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<RNG R, IntegrateMethod IM = IntegrateMethod::Simpson>
        [[nodiscard]] T compute_base(const Matrix auto& matA, Trv stepsize, int64_t numSample);

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] auto&& getGMRES(this auto&& self) noexcept { return self.gmres; }
        [[nodiscard]] size_t getOrder() const noexcept { return gmres.getLength(); }
    };

    template<Scalar T>
    SparseDet<T>::SparseDet(GMRES<T> gmres_)
            : gmres(std::move(gmres_)), precond(getOrder()), matD(getOrder()), rand(getOrder()), solX(getOrder()) {
        gmres.setTolerance(std::sqrt(std::numeric_limits<T>::epsilon()));
    }

    template<Scalar T>
    template<RNG R, IntegrateMethod IM>
    T SparseDet<T>::compute_base(const Matrix auto& matA, Trv stepsize, int64_t numSample) {
        assert(matA.getOrder() == getOrder());
        matD.diag() = matA.diag();

        Integrate<IM, Tr, 1> integral({{0}, {1}}, stepsize);
        T result = 0;
        for (int64_t i = 0; i < numSample; ++i) {
            rand.template random_rademacher<R>();
            const T sample = integral.solve([this, matA](Tr integralT) {
                precond.diag() = integralT + (Trv(1) - integralT) * matD.diag();
                for (auto& elem : precond.diag())
                    if (abs(elem) < Trv(1E-4)) // Select by experience
                        elem = 1;
                precond = precond.inv();

                auto m = integralT * matA + (Trv(1) - integralT) * matD;
                solX.zeros();
                gmres.solve(m * precond, rand, solX);
                return (rand.conjugate() * (matA - matD) * (precond * solX)).real();
            });
            result.toNextMean(i, sample);
        }
        return matD.lnAbsDet() + result;
    }

    template<Scalar T>
    void SparseDet<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        gmres.swap(obj.gmres);
        precond.swap(obj.precond);
        matD.swap(obj.matD);
        rand.swap(obj.rand);
        solX.swap(obj.solX);
    }
}
