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
#include <iostream>
#include <cassert>
#include <unistd.h>
#include <sstream>
#include <fstream>
#include <fcntl.h>
#include "Physica/Core/Interface/VaspWarpper.h"
#include "Physica/Core/Exception/BadFileFormatException.h"
#include "Physica/Utils/Unix/UnixHelper.h"

namespace Physica::Core {
    const char* VaspWarpper::errorMsg = "[Error]: Vasp finished with a non zero exit code";

    VaspWarpper::VaspWarpper()
            : vaspWorkingDir()
            , logFilePath()
            , core() {}

    VaspWarpper::VaspWarpper(size_t core_,
                             std::string pathToVasp_,
                             std::string workingDir,
                             std::string logFilePath_,
                             Poscar poscar_)
            : pathToVasp(std::move(pathToVasp_))
            , vaspWorkingDir(std::move(workingDir))
            , logFilePath(std::move(logFilePath_))
            , core(core_)
            , poscar(std::move(poscar_)) {
        std::ofstream fout((vaspWorkingDir + std::string("/POSCAR")).c_str(), std::ios_base::out | std::ios_base::trunc);
        fout << poscar;
        execute();
    }

    VaspWarpper& VaspWarpper::operator=(VaspWarpper vasp) noexcept {
        swap(vasp);
        return *this;
    }

    void VaspWarpper::execute() {
        future = Parallel::ProcessExecutor::schedule([=]() {
            int standardErr = dup(STDERR_FILENO);
            if (!logFilePath.empty()) {
                int log_fd = open(logFilePath.c_str()
                        , O_WRONLY | O_TRUNC | O_CREAT
                        , S_IRUSR | S_IWUSR);
                dup2(log_fd, STDOUT_FILENO);
                dup2(log_fd, STDERR_FILENO);
            }
            [[maybe_unused]] int err = chdir(vaspWorkingDir.c_str());
            /* Execute */ {
                constexpr size_t bufferLength = 21;
                char coreStr[bufferLength];
                [[maybe_unused]] int count = sprintf(coreStr, "%ld", core);
                assert(0 <= count && static_cast<size_t>(count) < bufferLength);
                execlp("mpirun", "mpirun", "-n", coreStr, pathToVasp.c_str(), nullptr);
            }
            dup2(standardErr, 2);
            perror("[Error]: Failed to execute VASP");
            _exit(EXIT_FAILURE);
        });
    }

    typename VaspWarpper::ScalarType VaspWarpper::getPress() const {
        future.wait(errorMsg);
        const char* tempfile = tmpnam(nullptr);
        const std::string command = std::string("grep 'in kB' ") +
                                    vaspWorkingDir +
                                    std::string("/OUTCAR | tr -s ' ' | cut -d ' ' -f 4,5,6 >") +
                                    std::string(tempfile);
        [[maybe_unused]] int err = system(command.c_str());
        ScalarType press_x, press_y, press_z;
        /* Read data */ {
            std::ifstream fin(tempfile);
            fin >> press_x >> press_y >> press_z;
        }
        unlink(tempfile);
        return (press_x + press_y + press_z) / 3.0f;
    }

    Outcar VaspWarpper::getOutcar() const {
        std::string path = vaspWorkingDir + std::string("/OUTCAR");
        return Outcar(path.c_str(), poscar.getAtomCount());
    }

    void VaspWarpper::swap(VaspWarpper& vasp) noexcept {
        pathToVasp.swap(vasp.pathToVasp);
        vaspWorkingDir.swap(vasp.vaspWorkingDir);
        logFilePath.swap(vasp.logFilePath);
        std::swap(core, vasp.core);
        poscar.swap(vasp.poscar);
        future.swap(vasp.future);
    }
}
