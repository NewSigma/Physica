/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Physics/ManyBody/ReprSpace/PauliMatrix.h"
#include "J1J2Matrix.h"

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof_tx<M, J1J2Matrix>)
    class GEMV<M, V> : public RValueVector<GEMV<M, V>> {
        using This = GEMV<M, V>;
        using Base = RValueVector<This>;
        using M1 = std::remove_cvref<M>::type;
        constexpr static int NumSite = M1::NumSite;
    protected:
        using typename Base::T;
        using typename Base::Trv;
    private:
        decay_rvalue_t<M> mat;
        decay_rvalue_t<V> vec;
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

        [[nodiscard]] T calc(size_t) const { noImpl(); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return mat.getRow(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isFastAssign() noexcept { return true; }
        [[nodiscard]] __host__ __device__ consteval static size_t getSizeAtCompile() noexcept { return M1::getRowAtCompile(); }
    private:
        /* Getters */
        [[nodiscard]] const auto& getRepr() const noexcept { return mat.getRepr(); }
        [[nodiscard]] auto getCouplingJ1() const noexcept { return mat.getCouplingJ1(); }
        [[nodiscard]] auto getCouplingJ2() const noexcept { return mat.getCouplingJ2(); }
    };

    template<Matrix M, Vector V> requires(instanceof_tx<M, J1J2Matrix>)
    GEMV<M, V>::GEMV(M&& mat_, V&& vec_) : mat(std::forward<M>(mat_)), vec(std::forward<V>(vec_)) {
        assert(mat.getCol() == vec.getLength());
    }

    template<Matrix M, Vector V> requires(instanceof_tx<M, J1J2Matrix>)
    template<ExecutePolicy P>
    void GEMV<M, V>::assign(Vector auto& target) const {
        target.assert_assign(*this);
        target.zeros();

        const auto& repr = getRepr();
        for (size_t i = 0; i < getLength(); ++i) {
            using enum PauliIndex;
            auto psiR = repr[i];
            /* Z component */ {
                int numAntiSpin = 0;
                int numNextAntiSpin = 0;
                for (int site = 0; site < NumSite; ++site) {
                    mat.forNeighSites([&, psiR](int from, int to) noexcept {
                        numAntiSpin += psiR[from] == psiR[to] ? 1 : -1;
                    }, site);

                    mat.forNNeighSites([&, psiR](int from, int to) noexcept {
                        numNextAntiSpin += psiR[from] == psiR[to] ? 1 : -1;
                    }, site);
                }
                target[i] += Trv(0.25) * vec.calc(i) * (getCouplingJ1() * T(numAntiSpin) + getCouplingJ2() * T(numNextAntiSpin));
            }
            /* XY component */ {
                for (int site = 0; site < NumSite; ++site) {
                    mat.forNeighSites([&, psiR](int from, int to) noexcept {
                        if (psiR[from] != psiR[to]) {
                            auto psiC = psiR.flip(from).flip(to);
                            target[repr[psiC]] += vec.calc(i) * getCouplingJ1() * T(0.5);
                        }
                    }, site);

                    mat.forNNeighSites([&, psiR](int from, int to) noexcept {
                        if (psiR[from] != psiR[to]) {
                            auto psiC = psiR.flip(from).flip(to);
                            target[repr[psiC]] += vec.calc(i) * getCouplingJ2() * T(0.5);
                        }
                    }, site);
                }
            }
        }
    }

    template<Matrix M, Vector V> requires(instanceof_tx<M, J1J2Matrix>)
    auto&& GEMV<M, V>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat);
    }

    template<Matrix M, Vector V> requires(instanceof_tx<M, J1J2Matrix>)
    auto&& GEMV<M, V>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), V>(self.vec);
    }
}

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof_tx<M, J1J2Matrix>)
    class Traits<GEMV<M, V>> {
        using T1 = std::remove_cvref_t<M>::ScalarType;
        using T2 = std::remove_cvref_t<V>::ScalarType;
    public:
        using ScalarType = Internal::BinaryScalarOpRtnTy<T1, T2>::Type;
    };
}
