/*
 * Copyright 2022-2026 Weibo He.
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

#include <fcntl.h>
#include "Physica/Core/Physics/SolidState/VASP/Outcar.h"
#include "Physica/Core/Parallel/Executor/ProcessExecutor.h"
#include "Physica/Core/Physics/MD/MDCell.h"
#include "Physica/Core/Utils/Unix/TempDir.h"

namespace Physica {
    template<Scalar T>
    class VASPModel {
        using This = VASPModel<T>;
        using MDCellType = MDCell<T>;
        using WorkingDirType = TempDir<15>;
    private:
        std::string pathToVasp;
        WorkingDirType workingDir;
        int logFd;
        Array<size_t> numOfEachType;
        unsigned int numMPIProcess;
    public:
        VASPModel() = default;
        VASPModel(
            std::string pathToVasp_,
            const char* pathToINCAR,
            const char* pathToPOTCAR,
            const char* pathToKpoints,
            Array<size_t> numOfEachType_,
            unsigned int numMPIProcess_);
        VASPModel(const This&) = delete;
        VASPModel(This&&) noexcept = default;
        ~VASPModel() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force(const MDCellType& cell) const;
        template<ExecutePolicy P>
        void forceAsync(const MDCellType& cell, Vector auto& result) const noexcept;
        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force_short(const MDCellType& cell) const { return force<P>(cell); }
        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force_long(const MDCellType& cell) const { return VectorND<T>(cell.getDOF(), 0); }
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const WorkingDirType& getWorkingDir() const noexcept { return workingDir; }
        [[nodiscard]] size_t getNumParticle() const noexcept;
        /* Static members */
        [[nodiscard]] consteval static bool isPeriodBoundary() noexcept { return true; }
    private:
        ProcessFuture run_vasp() const;
    };

    template<Scalar T>
    VASPModel<T>::VASPModel(std::string pathToVasp_,
                                     const char* pathToINCAR,
                                     const char* pathToPOTCAR,
                                     const char* pathToKPOINTS,
                                     Array<size_t> numOfEachType_,
                                     unsigned int numMPIProcess_)
            : pathToVasp(std::move(pathToVasp_))
            , workingDir("/tmp/tmpXXXXXX")
            , numOfEachType(std::move(numOfEachType_))
            , numMPIProcess(numMPIProcess_) {
        auto path = makePath("%s/log", workingDir.getName());
        logFd = open(path.data(), O_WRONLY | O_TRUNC | O_CREAT, S_IRUSR | S_IWUSR);
        if (logFd == -1)
            throw SystemException();

        path = makePath("%s/INCAR", workingDir.getName());
        copyFile(pathToINCAR, path.data());
        path = makePath("%s/POTCAR", workingDir.getName());
        copyFile(pathToPOTCAR, path.data());
        path = makePath("%s/KPOINTS", workingDir.getName());
        copyFile(pathToKPOINTS, path.data());
    }

    template<Scalar T>
    template<ExecutePolicy P>
    VectorND<T> VASPModel<T>::force(const MDCellType& cell) const {
        VectorND<T> result(cell.getNumParticle());
        forceAsync<P>(cell, result);
        return result;
    }

    template<Scalar T>
    template<ExecutePolicy P>
    void VASPModel<T>::forceAsync(const MDCellType& cell, Vector auto& result) const noexcept {
        try {
            /* Make POSCAR */ {
                const auto path = makePath("%s/POSCAR", workingDir.getName());
                std::ofstream os(path.data(), std::ios_base::out | std::ios_base::trunc);
                os << '\n';
                os << PhyConst<AU>::bohrToAngstorm(1) << '\n';
                os << cell.getLattice();
                for (size_t elem : numOfEachType)
                    os << elem << ' ';
                os << "\nCartesian\n";
                os << cell.getPos();
            }
            const int error = run_vasp().wait();
            if (error != 0) [[unlikely]]
                throw std::runtime_error("[Error]: VASP finished with non zero exit code");
            const auto path = makePath("%s/OUTCAR", workingDir.getName());
            const Outcar outcar(path.data(), getNumParticle());
            result.getDerived() = outcar.getForce();
        }
        catch (std::exception& e) {
            fprintf(stderr, "%s\n", e.what());
            exit(EXIT_FAILURE);
        }
    }

    template<Scalar T>
    void VASPModel<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        pathToVasp.swap(obj.pathToVasp);
        workingDir.swap(obj.workingDir);
        std::swap(logFd, obj.logFd);
        numOfEachType.swap(obj.numOfEachType);
        std::swap(numMPIProcess, obj.numMPIProcess);
    }

    template<Scalar T>
    size_t VASPModel<T>::getNumParticle() const noexcept {
        size_t result = 0;
        for (size_t elem : numOfEachType)
            result += elem;
        return result;
    }

    template<Scalar T>
    ProcessFuture VASPModel<T>::run_vasp() const {
        return ProcessExecutor::schedule([&]() {
            const int standardErr = dup(STDERR_FILENO);
            dup2(logFd, STDOUT_FILENO);
            dup2(logFd, STDERR_FILENO);
            [[maybe_unused]] int err = chdir(workingDir.getName());
            /* Execute */ {
                constexpr int BufferSize = 16; // 16 is enough for unsigned int
                std::array<char, BufferSize> numProcess{};
                [[maybe_unused]] const auto count = std::format_to_n(numProcess.data(), numProcess.size(), "%d", numMPIProcess).size;
                assert(0 <= count && count < BufferSize && "[Error]: Unexpected bad printf");
                execlp("mpirun", "mpirun", "-n", numProcess, pathToVasp.c_str(), nullptr);
            }
            dup2(standardErr, STDERR_FILENO);
            perror("[Error]: Failed to execute VASP");
            _exit(EXIT_FAILURE);
        });
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<VASPModel<T>> {
    public:
        constexpr static bool IsContractable = false;
    };
}
