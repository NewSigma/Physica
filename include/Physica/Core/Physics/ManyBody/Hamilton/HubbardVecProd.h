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
#pragma once

#include "HubbardMatrix.h"

namespace Physica {
    template<Scalar T0, Representation U, BoundaryCond BC, Vector V>
    class HubbardVecProd : public RValueVector<HubbardVecProd<T0, U, BC, V>> {
        using MatrixType = HubbardMatrix<T0, U, BC>;
        using This = HubbardVecProd<T0, U, BC, V>;
        using Base = RValueVector<This>;
        using FFT1D = MatrixType::FFT1D;
        using StateType = U::StateType;
        constexpr static unsigned int Dim = MatrixType::Dim;
        constexpr static unsigned int NumSite = MatrixType::NumSite;
        constexpr static unsigned int SiteDOF = MatrixType::SiteDOF;
        constexpr static bool IsTransInvariant = MatrixType::IsTransInvariant;
    public:
        using typename Base::ScalarType;
    protected:
        using typename Base::T;
        using typename Base::Tr;
        using typename Base::Tv;
    private:
        const MatrixType& mat;
        const V& vec;
    public:
        HubbardVecProd(const MatrixType& mat_, const V& vec_);
        HubbardVecProd(const This&) = default;
        HubbardVecProd(This&&) noexcept = default;
        ~HubbardVecProd() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& target) const;

        [[nodiscard]] T calc(size_t index) const;
        [[nodiscard]] Tv calc_value(size_t index) const;
        /* Getters */
        [[nodiscard]] size_t getLength() const { return mat.getRow(); }
        [[nodiscard]] const MatrixType& getLHS() const noexcept { return mat; }
        [[nodiscard]] const V& getRHS() const noexcept { return vec; }
    private:
        void sumHopping(Vector auto& target, FFT1D& fft, T factor, StateType psi) const;
        void dotImpl(Vector auto& target, T factor, size_t index) const;
        /* Getters */
        [[nodiscard]] T getHoppingT() const noexcept { return mat.getHoppingT(); }
        [[nodiscard]] T getRepelU() const noexcept { return mat.getRepelU(); }
        [[nodiscard]] const U& getRepr() const noexcept { return mat.getRepr(); }
    };

    template<Scalar T0, Representation U, BoundaryCond BC, Vector V>
    HubbardVecProd<T0, U, BC, V>::HubbardVecProd(const MatrixType& mat_, const V& vec_) : mat(mat_), vec(vec_) {
        assert(mat.getCol() == vec.getLength());
    }

    template<Scalar T0, Representation U, BoundaryCond BC, Vector V>
    template<ExecutePolicy P>
    void HubbardVecProd<T0, U, BC, V>::assign(Vector auto& target) const {
        assert(target.getLength() == getLength() && "[Error]: Dimensions do not match");
        if constexpr (P == Thread) {
            parallel_for<Thread>([&](size_t i) {
                target[i] = calc(i);
            }, getLength(), 0).wait();
        }
        else {
            target = Tr(0);
            const size_t length = getLength();
            for (size_t i = 0; i < length; ++i)
                dotImpl(target, vec.calc(i), i);
        }
    }

    template<Scalar T0, Representation U, BoundaryCond BC, Vector V>
    auto HubbardVecProd<T0, U, BC, V>::calc(size_t index) const -> T {
        static_assert(!IsTransInvariant && "[Error]: Not implemented");
        const T hop = -mat.getHoppingT();
        const auto state = getRepr()[index];
        T result = 0;
        int numRepel = 0;
        for (int site = 0; site < int(NumSite); ++site) {
            const bool upOccupy1 = state.isUpOccupy(site);
            const bool downOccupy1 = state.isDownOccupy(site);
            numRepel += upOccupy1 && downOccupy1;

            mat.forNeighSites([this, &result, &state, hop, upOccupy1, downOccupy1](int site, int site1) {
                const bool upOccupy2 = state.isUpOccupy(site1);
                const bool downOccupy2 = state.isDownOccupy(site1);
                const Vector2D<T> phases = mat.calcBoundaryPhase(site, site1);
                if (upOccupy1 != upOccupy2) {
                    const T hopUp = hop * Tr(state.hopUpSign(site, site1));
                    const size_t index1 = getRepr()[upOccupy1 ? state.hopUp(site, site1) : state.hopUp(site1, site)];
                    result += vec.calc(index1) * phases[!upOccupy1] * (upOccupy1 ? hopUp : -hopUp);
                }

                if (downOccupy1 != downOccupy2) {
                    const T hopDown = hop * Tr(state.hopDownSign(site, site1));
                    const size_t index1 = getRepr()[downOccupy1 ? state.hopDown(site, site1) : state.hopDown(site1, site)];
                    result += vec.calc(index1) * phases[!downOccupy1] * (downOccupy1 ? hopDown : -hopDown);
                }
            }, site);
        }
        result += vec.calc(index) * (mat.getRepelU() * Tr(numRepel));
        return result;
    }

    template<Scalar T0, Representation U, BoundaryCond BC, Vector V>
    auto HubbardVecProd<T0, U, BC, V>::calc_value(size_t index) const -> Tv {
        return calc(index).value();
    }

    template<Scalar T0, Representation U, BoundaryCond BC, Vector V>
    void HubbardVecProd<T0, U, BC, V>::sumHopping(Vector auto& target, FFT1D& fft, T factor, StateType psi) const {
        if (psi.isVacuum())
            return;
        const auto reducedPsi = psi.transReduce();
        auto& rSpace = fft.getRSpace();
        int sign = 1;
        for (int i = 0; i < int(NumSite); ++i) {
            rSpace[i] = Tr(reducedPsi == psi ? sign : 0);
            sign *= psi.lShiftSign();
            psi <<= 1;
        }
        FFT1D::transform(mat.planProvider, fft);
        const auto& repr = getRepr();
        const size_t index = repr[reducedPsi];
        target[index] += fft.getKSpace()[repr.getReducedK()] * sqrt(Tr(repr.getPeriods()[index])) * factor;
    }

    template<Scalar T0, Representation U, BoundaryCond BC, Vector V>
    void HubbardVecProd<T0, U, BC, V>::dotImpl(Vector auto& target, T factor, size_t index) const {
        const auto state = getRepr()[index];
        /* On site contribution */ {
            int numRepel = 0;
            for (int site = 0; site < int(NumSite); ++site)
                numRepel += state.isUpOccupy(site) && state.isDownOccupy(site);
            target[index] += factor * (getRepelU() * Tr(numRepel));
        }

        if constexpr (IsTransInvariant) {
            static_assert(Dim == 1 && "[Error]: Not implemented");
            static_assert(BC == BoundaryCond::PBC, "[Error]: Not implemented");
            const Tr normalizer = sqrt(Tr(getRepr().getPeriods()[index])) / Tr(NumSite);
            const T hop = -factor * normalizer * getHoppingT();

            auto fft = FFT1D::makeEmptyFFT(NumSite);
            for (int site = 0; site < int(NumSite); ++site) {
                const auto site1 = (site + 1) % NumSite;
                const T hopUp = hop * Tr(state.hopUpSign(site, site1));
                const T hopDown = hop * Tr(state.hopDownSign(site, site1));
                sumHopping(target, fft, hopUp, state.hopUp(site, site1));
                sumHopping(target, fft, -hopUp, state.hopUp(site1, site));
                sumHopping(target, fft, hopDown, state.hopDown(site, site1));
                sumHopping(target, fft, -hopDown, state.hopDown(site1, site));
            }
        }
        else {
            const T hop = -factor * getHoppingT();
            for (int site = 0; site < int(NumSite); ++site) {
                const bool upOccupy1 = state.isUpOccupy(site);
                const bool downOccupy1 = state.isDownOccupy(site);
                mat.forNeighSites([this, &target, &state, hop, upOccupy1, downOccupy1](int site, int site1) noexcept {
                    const auto& repr = getRepr();
                    const bool upOccupy2 = state.isUpOccupy(site1);
                    const bool downOccupy2 = state.isDownOccupy(site1);
                    const Vector2D<T> phases = mat.calcBoundaryPhase(site, site1);
                    if (upOccupy1 != upOccupy2) {
                        const size_t index = repr[upOccupy1 ? state.hopUp(site, site1) : state.hopUp(site1, site)];
                        const bool sign = state.hopUpSign(site, site1) == 1;
                        target[index] += ((upOccupy1 == sign) ? hop : -hop) * phases[!upOccupy1];
                    }

                    if (downOccupy1 != downOccupy2) {
                        const size_t index = repr[downOccupy1 ? state.hopDown(site, site1) : state.hopDown(site1, site)];
                        const bool sign = state.hopDownSign(site, site1) == 1;
                        target[index] += (downOccupy1 == sign) ? hop : -hop * phases[!downOccupy1];
                    }
                }, site);
            }
        }
    }
}

namespace Physica {
    template<Scalar T0, Representation U, BoundaryCond BC, Vector V>
    class Traits<HubbardVecProd<T0, U, BC, V>> {
        using MatrixType = HubbardMatrix<T0, U, BC>;
        using T1 = V::ScalarType;
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<T0, T1>::Type;
        constexpr static size_t SizeAtCompile = MatrixType::RowAtCompile;
        constexpr static bool FastAssign = true;
        constexpr static bool FastPacket = false;
    };
}
