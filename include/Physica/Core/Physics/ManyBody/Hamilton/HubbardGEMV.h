/*
 * Copyright 2024-2026 Weibo He.
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
    template<Matrix M, Vector V> requires(instanceof_ttx<M, HubbardMatrix>)
    class GEMV<M, V> : public RValueVector<GEMV<M, V>> {
        using This = GEMV<M, V>;
        using Base = RValueVector<This>;

        using M1 = std::remove_cvref<M>::type;
        using FFT1D = M1::FFT1D;
        using StateType = M1::StateType;
        constexpr static unsigned int Dim = M1::Dim;
        constexpr static unsigned int NumSite = M1::NumSite;
        constexpr static unsigned int SiteDOF = M1::SiteDOF;
        constexpr static bool IsTransInvariant = M1::IsTransInvariant;
    public:
        using typename Base::ScalarType;
    protected:
        using typename Base::T;
        using typename Base::Tr;
        using typename Base::Tv;
    private:
        LazyDestroy<M> mat;
        LazyDestroy<V> vec;
    public:
        GEMV(M&& mat_, V&& vec_);
        GEMV(const This&) = default;
        GEMV(This&&) noexcept = default;
        ~GEMV() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Vector auto& target) const;

        [[nodiscard]] T calc(size_t index) const;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isFastAssign() noexcept { return true; }
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return Dynamic; }
    private:
        void sumHopping(Vector auto& target, FFT1D& fft, T factor, StateType psi) const;
        void dotImpl(Vector auto& target, T factor, size_t index) const;
        /* Getters */
        [[nodiscard]] T getHoppingT() const noexcept { return mat.getHoppingT(); }
        [[nodiscard]] T getRepelU() const noexcept { return mat.getRepelU(); }
        [[nodiscard]] const auto& getRepr() const noexcept { return mat.getRepr(); }
    };

    template<Matrix M, Vector V> requires(instanceof_ttx<M, HubbardMatrix>)
    GEMV<M, V>::GEMV(M&& mat_, V&& vec_) : mat(std::forward<M>(mat_)), vec(std::forward<V>(vec_)) {
        assert(mat.getCol() == vec.getLength());
    }

    template<Matrix M, Vector V> requires(instanceof_ttx<M, HubbardMatrix>)
    template<ExecutePolicy P>
    void GEMV<M, V>::assign(Vector auto& target) const {
        target.assert_assign(*this);
        if constexpr (P == Thread) {
            parallel_for<Thread>([&](size_t i) {
                target[i] = calc(i);
            }, getLength(), 0).wait();
        }
        else {
            target.zeros();
            const size_t length = getLength();
            for (size_t i = 0; i < length; ++i)
                dotImpl(target, vec.calc(i), i);
        }
    }

    template<Matrix M, Vector V> requires(instanceof_ttx<M, HubbardMatrix>)
    auto GEMV<M, V>::calc(size_t index) const -> T {
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
                    result += vec.calc(index1) * phases[upOccupy1] * (upOccupy1 ? hopUp : -hopUp);
                }

                if (downOccupy1 != downOccupy2) {
                    const T hopDown = hop * Tr(state.hopDownSign(site, site1));
                    const size_t index1 = getRepr()[downOccupy1 ? state.hopDown(site, site1) : state.hopDown(site1, site)];
                    result += vec.calc(index1) * phases[downOccupy1] * (downOccupy1 ? hopDown : -hopDown);
                }
            }, site);
        }
        result += vec.calc(index) * (mat.getRepelU() * Tr(numRepel));
        return result;
    }

    template<Matrix M, Vector V> requires(instanceof_ttx<M, HubbardMatrix>)
    auto&& GEMV<M, V>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }

    template<Matrix M, Vector V> requires(instanceof_ttx<M, HubbardMatrix>)
    auto&& GEMV<M, V>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.vec);
    }

    template<Matrix M, Vector V> requires(instanceof_ttx<M, HubbardMatrix>)
    void GEMV<M, V>::sumHopping(Vector auto& target, FFT1D& fft, T factor, StateType psi) const {
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
        fft.transform(mat.getPlan(), fft);
        const auto& repr = getRepr();
        const size_t index = repr[reducedPsi];
        target[index] += fft.getKSpace()[repr.getReducedK()] * sqrt(Tr(repr.getPeriods()[index])) * factor;
    }

    template<Matrix M, Vector V> requires(instanceof_ttx<M, HubbardMatrix>)
    void GEMV<M, V>::dotImpl(Vector auto& target, T factor, size_t index) const {
        const auto state = getRepr()[index];
        /* On site contribution */ {
            int numRepel = 0;
            for (int site = 0; site < int(NumSite); ++site)
                numRepel += state.isUpOccupy(site) && state.isDownOccupy(site);
            target[index] += factor * (getRepelU() * Tr(numRepel));
        }

        if constexpr (IsTransInvariant) {
            static_assert(Dim == 1 && "[Error]: Not implemented");
            static_assert(Traits<M>::Boundary == BoundaryCond::PBC, "[Error]: Not implemented");
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
                        const bool sign = upOccupy1 == (state.hopUpSign(site, site1) == 1);
                        target[index] += (sign ? hop : -hop) * phases[upOccupy2];
                    }

                    if (downOccupy1 != downOccupy2) {
                        const size_t index = repr[downOccupy1 ? state.hopDown(site, site1) : state.hopDown(site1, site)];
                        const bool sign = downOccupy1 == (state.hopDownSign(site, site1) == 1);
                        target[index] += (sign ? hop : -hop) * phases[downOccupy2];
                    }
                }, site);
            }
        }
    }
}

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof_ttx<M, HubbardMatrix>)
    class Traits<GEMV<M, V>> {
        using M1 = std::remove_cvref<M>::type;
        using V1 = std::remove_cvref<V>::type;
        using T1 = M1::ScalarType;
        using T2 = V1::ScalarType;
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<T1, T2>::Type;
    };
}
