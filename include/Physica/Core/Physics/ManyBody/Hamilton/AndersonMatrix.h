/*
 * Copyright 2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomp/Tridiagonalization.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Physics/ManyBody/Hamilton/HamiltonMatrix.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/FermiRepr.h"

namespace Physica {
    /**
     * Reference:
     * [1] Gubernatis J, Kawashima N, Werner P. Quantum Monte Carlo Methods: Algorithms for Lattice Models. Cambridge University Press; 2016:214-217  
     */
    template<Scalar T, int NumSite, bool UseInversionSymm>
    class AndersonMatrix : public HamiltonMatrix<AndersonMatrix<T, NumSite, UseInversionSymm>> {
        using This = AndersonMatrix<T, NumSite, UseInversionSymm>;
        using Base = HamiltonMatrix<This>;
        using Repr = Traits<This>::ReprType;

        VectorND<T> conductE;
        VectorND<T> hybridV;
        T repelU;
        T chemMu;
        Repr repr;

        MatrixND<T> chainRepr;
    public:
        AndersonMatrix(T repelU, T chemMu, VectorND<T> conductE, VectorND<T> hybridV, Repr repr);
        AndersonMatrix(const This&) = default;
        AndersonMatrix(This&&) noexcept = default;
        ~AndersonMatrix() = default;
        /* Operators */
        This& operator=(This obj) { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] T calc(size_t row, size_t col) const;

        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] constexpr static int getNumSite() noexcept { return NumSite; }
        [[nodiscard]] T getRepelU() const noexcept { return repelU; }
        [[nodiscard]] T getChemMu() const noexcept { return chemMu; }
        [[nodiscard]] const auto& getRepr() const noexcept { return repr; }
        [[nodiscard]] T getOnSiteE(int site) const noexcept;
        [[nodiscard]] T getHoppingT(int bond) const noexcept;
    private:
        MatrixND<T> makeSingleParticleMatrix() const;
    };

    template<Scalar T, int NumSite, bool UseInversionSymm>
    AndersonMatrix<T, NumSite, UseInversionSymm>::AndersonMatrix(T repelU, T chemMu, VectorND<T> conductE, VectorND<T> hybridV, Repr repr)
            : conductE(std::move(conductE))
            , hybridV(std::move(hybridV))
            , repelU(repelU)
            , chemMu(chemMu)
            , repr(std::move(repr)) {
        assert((this->conductE.getLength() == this->hybridV.getLength()) && "[Error]: Conduct site number is not consistent");
        assert((this->conductE.getLength() == NumSite - 1) && "[Error]: NumSite it not consistent");
        Tridiagonalization<T> tri(makeSingleParticleMatrix());
        chainRepr = tri.getMatrixT();
    }

    template<Scalar T, int NumSite, bool UseInversionSymm>
    T AndersonMatrix<T, NumSite, UseInversionSymm>::calc(size_t row, size_t col) const {
        const auto psiR = repr[row];
        if (row == col) {
            T result = 0;
            for (int site = 0; site < NumSite; ++site) {
                bool occU = psiR.isUpOccupy(site);
                bool occD = psiR.isDownOccupy(site);
                if (site == 0)
                    result += (occU && occD) ? getRepelU() : T(0);
                result += T(int(occU) + int(occD)) * getOnSiteE(site);
            }
            return result;
        }

        const auto psiC = repr[col];
        T result = 0;
        for (int bond = 0; bond < NumSite - 1; ++bond) {
            int bond1 = bond + 1;
            const int signUp = psiC.hopUpSign(bond, bond1);
            const int signDown = psiC.hopDownSign(bond, bond1);
            int n1 = (psiR == psiC.hopUp(bond, bond1))
                   - (psiR == psiC.hopUp(bond1, bond));
            int n2 = (psiR == psiC.hopDown(bond, bond1))
                   - (psiR == psiC.hopDown(bond1, bond));
            result += T(n1 * signUp + n2 * signDown) * getHoppingT(bond);
        }
        return result;
    }

    template<Scalar T, int NumSite, bool UseInversionSymm>
    void AndersonMatrix<T, NumSite, UseInversionSymm>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        conductE.swap(obj.conductE);
        hybridV.swap(obj.hybridV);
        repelU.swap(obj.repelU);
        chemMu.swap(obj.chemMu);
        repr.swap(obj.repr);

        chainRepr.swap(obj.chainRepr);
    }

    template<Scalar T, int NumSite, bool UseInversionSymm>
    T AndersonMatrix<T, NumSite, UseInversionSymm>::getOnSiteE(int site) const noexcept {
        return chainRepr(site, site);
    }

    template<Scalar T, int NumSite, bool UseInversionSymm>
    T AndersonMatrix<T, NumSite, UseInversionSymm>::getHoppingT(int bond) const noexcept {
        assert(bond + 1 < getNumSite());
        return chainRepr(bond + 1, bond);
    }

    template<Scalar T, int NumSite, bool UseInversionSymm>
    MatrixND<T> AndersonMatrix<T, NumSite, UseInversionSymm>::makeSingleParticleMatrix() const {
        auto result = MatrixND<T>(NumSite);
        result.zeros();

        auto diag = result.diag();
        diag[0] = -chemMu;
        diag.tail(1) = conductE;
        result.col(0).tail(1) = hybridV;
        result.row(0).tail(1) = hybridV;
        return result;
    }
}

namespace Physica {
    template<Scalar T, int I, bool B>
    class Traits<AndersonMatrix<T, I, B>> : public Traits<HamiltonMatrix<AndersonMatrix<T, I, B>>> {
    public:
        using ScalarType = T;
        constexpr static int NumSite = I;
        constexpr static bool UseInversionSymm = B;
        constexpr static auto Boundary = BoundaryCond::OBC;

        using ReprType = FermiRepr<1, NumSite, UseInversionSymm>;
    };
}

#include "AndersonGEMV.h"
