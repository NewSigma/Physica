/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2024 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <fstream>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica {
    class PHYSICA_API PWscfOut {
        using This = PWscfOut;
        using ScalarType = float64;
        constexpr static size_t DefaultBufferSize = 1024; //1024 shall be enough

        std::ifstream fin;
        VectorND<ScalarType> force;
        Array<char> buffer;
    public:
        PWscfOut(const char* path, size_t numAtom);
        PWscfOut(const This&) = delete;
        PWscfOut(This&&) noexcept = default;
        ~PWscfOut() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        [[nodiscard]] VectorND<ScalarType> makeTotalForces();
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumAtom() const noexcept { return force.getLength() / 3; }
        [[nodiscard]] const VectorND<ScalarType>& getForce() const noexcept { return force; }
    private:
        void readForce();
    };
}
