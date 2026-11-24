/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2024 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <fstream>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica {
    class PHYSICA_API Outcar {
        using This = Outcar;
        using ScalarType = float64;

        VectorND<ScalarType> force;
        ScalarType internalEnergy;
    public:
        Outcar(const char* path, unsigned int numAtom);
        Outcar(const This&) = default;
        Outcar(This&&) noexcept = default;
        ~Outcar() = default;
        /* Operators */
        Outcar& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumAtom() const noexcept { return force.getLength() / 3; }
        [[nodiscard]] const VectorND<ScalarType>& getForce() const noexcept { return force; }
        [[nodiscard]] ScalarType getInternalEnergy() const noexcept { return internalEnergy; }
    private:
        void readForce(std::ifstream& fin, Array<char>& buffer);
        void readEnergy(std::ifstream& fin, Array<char>& buffer);
    };
}
