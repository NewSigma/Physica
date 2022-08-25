/*
 * Copyright 2022 WeiBo He.
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

#include <cstdlib>
#include <memory>
#include "Physica/Core/Parallel/SubProcess.h"
#include "Physica/Core/IO/Poscar.h"

namespace Physica::Core {
    class VaspWarpper final : public Parallel::SubProcess {
        using Base = Parallel::SubProcess;
        using ScalarType = Scalar<Float, false>;
    private:
        std::string pathToVasp;
        std::string vaspWorkingDir;
        std::string logFilePath;
        size_t core;
        Poscar poscar;
    public:
        VaspWarpper();
        VaspWarpper(size_t core_, std::string pathToVasp_, std::string workingDir, std::string logFilePath_, Poscar poscar_);
        VaspWarpper(const VaspWarpper&) = delete;
        VaspWarpper(VaspWarpper&& vasp) noexcept;
        ~VaspWarpper() = default;
        /* Operators */
        VaspWarpper& operator=(VaspWarpper vasp) noexcept;
        /* Getters */
        [[nodiscard]] const std::string& getWorkingDir() const noexcept { return vaspWorkingDir; }
        [[nodiscard]] ScalarType getEnergy() const;
        [[nodiscard]] ScalarType getPress() const;
        [[nodiscard]] Vector<ScalarType> getForce() const;
        /* Helpers */
        void swap(VaspWarpper& vasp) noexcept;
    private:
        void run() const;
    };
}
