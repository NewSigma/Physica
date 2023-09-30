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
#include "Physica/Core/IO/Poscar.h"
#include "Physica/Core/IO/Outcar.h"
#include "Physica/Core/Parallel/Executor/ProcessExecutor.h"

namespace Physica::Core {
    class VaspWarpper final {
        using ScalarType = Scalar<Float>;
        static const char* errorMsg;
    private:
        std::string pathToVasp;
        std::string vaspWorkingDir;
        std::string logFilePath;
        size_t core;
        Poscar<ScalarType> poscar;
        mutable ProcessFuture future;
    public:
        VaspWarpper();
        VaspWarpper(size_t core_, std::string pathToVasp_, std::string workingDir, std::string logFilePath_, Poscar<ScalarType> poscar_);
        VaspWarpper(const VaspWarpper&) = delete;
        VaspWarpper(VaspWarpper&&) noexcept = default;
        ~VaspWarpper() = default;
        /* Operators */
        VaspWarpper& operator=(VaspWarpper vasp) noexcept;
        /* Operations */
        void execute();
        /* Getters */
        [[nodiscard]] const std::string& getWorkingDir() const noexcept { return vaspWorkingDir; }
        [[nodiscard]] ScalarType getPress() const;
        [[nodiscard]] Outcar getOutcar() const;
        /* Helpers */
        void swap(VaspWarpper& vasp) noexcept;
    private:
        friend class Test;
    };
}
