/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <print>
#include "Physica/Core/Scalar/Real.h"
#include "Physica/Core/Parallel/MPI.h"

using namespace Physica;

namespace {
    template<FloatPrec Prec>
    Real<Prec> stripSignificandImpl(Real<Prec> x) noexcept {
        assert(!x.isSubNormal() && "[Error]: Input must be a normal positive value");
        int exp = 0;
        std::frexp(x.toMachine(), &exp);
        return Real<Prec>(std::ldexp(Real<Prec>(1).toMachine(), exp - 1));
    }
}

auto Real<Float32>::stripSignificand() const noexcept -> Real {
    return stripSignificandImpl(*this);
}

void Real<Float32>::dump() const noexcept {
    std::println("{}", *this);
}

auto Real<Float64>::stripSignificand() const noexcept -> Real {
    return stripSignificandImpl(*this);
}

void Real<Float64>::dump() const noexcept {
    std::println("{}", *this);
}

auto Real<Float32>::dtype_mpi() noexcept -> Handle<HandleType::MPI_Dtype> {
    return MPI::dtype<float>();
}

auto Real<Float64>::dtype_mpi() noexcept -> Handle<HandleType::MPI_Dtype> {
    return MPI::dtype<double>();
}
