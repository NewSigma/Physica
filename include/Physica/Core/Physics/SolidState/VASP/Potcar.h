/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <fstream>
#include "Physica/Core/Scalar/Real.h"

namespace Physica {
    class PHYSICA_API Potcar {
        using ScalarType = float32;
        ScalarType mass;
        unsigned char numValenceElectron;
    public:
        explicit Potcar(const char* path);
        Potcar(const Potcar&) = default;
        Potcar(Potcar&&) noexcept = default;
        ~Potcar() = default;
        /* Operators */
        Potcar& operator=(Potcar obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(Potcar& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] ScalarType getMass() const noexcept { return mass; }
        [[nodiscard]] unsigned int getNumValenceElectron() const noexcept { return numValenceElectron; }
    private:
        void readFile(std::ifstream& fin);
    };
}
