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
#include <cassert>
#include <unistd.h>
#include <sstream>
#include <fstream>
#include <fcntl.h>
#include "Physica/Core/Interface/VaspWarpper.h"
#include "Physica/Utils/Unix/UnixHelper.h"

namespace Physica::Core {
    VaspWarpper::VaspWarpper()
            : Base([this]() { this->run(); })
            , vaspWorkingDir()
            , logFilePath()
            , core() {}

    VaspWarpper::VaspWarpper(size_t core_,
                             std::string pathToVasp_,
                             std::string workingDir,
                             std::string logFilePath_,
                             const Poscar& poscar)
            : Base([this]() { this->run(); })
            , pathToVasp(std::move(pathToVasp_))
            , vaspWorkingDir(std::move(workingDir))
            , logFilePath(std::move(logFilePath_))
            , core(core_) {
        std::ofstream fout((std::string(workingDir) + std::string("/POSCAR")).c_str(), std::ios_base::out | std::ios_base::trunc);
        fout << poscar;
    }
    /**
     * \param poscarId
     * The id of poscar to be calculated.
     * \param groupNum
     * The number of current group. Equals to groupId + 1.
     */
    VaspWarpper::VaspWarpper(size_t core_, size_t poscarId)
            : SubProcess([this]() { this->run(); })
            , vaspWorkingDir()
            , logFilePath()
            , core(core_) {
        std::ostringstream ostr{};
        ostr << "runvasp/" << poscarId;
        vaspWorkingDir = ostr.str();
        ostr.str("");
        ostr << "tmpdata/output_" << poscarId;
        logFilePath = ostr.str();
        /* Move POSCAR */ {
            auto fromPoscar = Utils::makePath("POSCAR_%ld", poscarId);
            auto toPoscar = Utils::makePath("runvasp/%ld/POSCAR", poscarId);
            rename(fromPoscar.get(), toPoscar.get());
        }
    }

    VaspWarpper::VaspWarpper(VaspWarpper&& vasp) noexcept
            : SubProcess(std::move(vasp))
            , pathToVasp(std::move(vasp.pathToVasp))
            , vaspWorkingDir(std::move(vasp.vaspWorkingDir))
            , logFilePath(std::move(vasp.logFilePath))
            , core(vasp.core) {}

    VaspWarpper& VaspWarpper::operator=(VaspWarpper vasp) noexcept {
        swap(vasp);
        return *this;
    }

    void VaspWarpper::run() const {
        int standardErr = dup(2);
        if (!logFilePath.empty()) {
            int log_fd = open(logFilePath.c_str()
                    , O_WRONLY | O_TRUNC | O_CREAT
                    , S_IRUSR | S_IWUSR);
            dup2(log_fd, 1);
            dup2(log_fd, 2);
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
    }

    float VaspWarpper::getEnergy() const {
        const char* tempfile = tmpnam(nullptr);
        const std::string command = std::string("grep energy ") +
                                    vaspWorkingDir +
                                    std::string("/OUTCAR | tail -n 1 | tr -s ' ' | cut -d ' ' -f 5 >") +
                                    std::string(tempfile);
        [[maybe_unused]] int err = system(command.c_str());
        float result;
        /* Read data */ {
            std::ifstream fin(tempfile);
            fin >> result;
        }
        unlink(tempfile);
        return result;
    }

    float VaspWarpper::getPress() const {
        const char* tempfile = tmpnam(nullptr);
        const std::string command = std::string("grep 'in kB' ") +
                                    vaspWorkingDir +
                                    std::string("/OUTCAR | tr -s ' ' | cut -d ' ' -f 4,5,6 >") +
                                    std::string(tempfile);
        [[maybe_unused]] int err = system(command.c_str());
        float press_x, press_y, press_z;
        /* Read data */ {
            std::ifstream fin(tempfile);
            fin >> press_x >> press_y >> press_z;
        }
        unlink(tempfile);
        return (press_x + press_y + press_z) / 3.0f;
    }

    void VaspWarpper::swap(VaspWarpper& vasp) noexcept {
        pathToVasp.swap(vasp.pathToVasp);
        vaspWorkingDir.swap(vasp.vaspWorkingDir);
        logFilePath.swap(vasp.logFilePath);
        std::swap(core, vasp.core);
    }
}
