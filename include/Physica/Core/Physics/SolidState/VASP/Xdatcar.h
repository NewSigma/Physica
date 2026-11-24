/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <fstream>
#include <utility>
#include "Poscar.h"

namespace Physica {
    class PHYSICA_API Xdatcar {
        using ScalarType = float64;
    private:
        Poscar<ScalarType> data;
        std::ifstream fin;
        uint64_t stepNum;
        bool init;
    public:
        Xdatcar(std::ifstream fin_);
        Xdatcar(const Xdatcar&) = delete;
        Xdatcar(Xdatcar&&) noexcept = default;
        ~Xdatcar() = default;
        /* Operators */
        Xdatcar& operator=(Xdatcar xdatcar) noexcept;
        /* Operations */
        bool step();
        void swap(Xdatcar& __restrict xdatcar) noexcept;
        /* Getters */
        [[nodiscard]] const Poscar<ScalarType>& getCurrent() const noexcept { return data; }
        [[nodiscard]] uint64_t getStep() const noexcept { return stepNum; }
    };
}
