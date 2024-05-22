/*
 * Copyright 2022-2024 WeiBo He.
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
#include "Physica/Core/IO/VASP/Poscar.h"
#include "Physica/Core/IO/Outcar.h"
#include "Physica/Core/Parallel/Executor/ProcessExecutor.h"
#include "Physica/Core/Physics/MD/MDCell.h"
#include "Physica/Utils/Unix/TempDir.h"
#include "Physica/Utils/Unix/UnixHelper.h"

namespace Physica::Core {
    template<class ScalarType> class VASPModel;

    namespace Internal {
        template<class ScalarType>
        class Traits<VASPModel<ScalarType>> {
        public:
            constexpr static bool IsPeriodBoundary = true;
            constexpr static bool IsContractable = false;
        };
    }

    template<class ScalarType>
    class VASPModel {
        using MDCellType = MDCell<ScalarType>;
        using WorkingDirType = Utils::TempDir<15>;
    private:
        std::string pathToVasp;
        WorkingDirType workingDir;
        int logFd;
        Utils::Array<size_t> numOfEachType;
        unsigned int numMPIProcess;
    public:
        VASPModel() = default;
        VASPModel(
            std::string pathToVasp_,
            const char* pathToINCAR,
            const char* pathToPOTCAR,
            const char* pathToKpoints,
            Utils::Array<size_t> numOfEachType_,
            unsigned int numMPIProcess_);
        VASPModel(const VASPModel&) = delete;
        VASPModel(VASPModel&&) noexcept = default;
        ~VASPModel() = default;
        /* Operators */
        VASPModel& operator=(VASPModel obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force(const MDCellType& cell) const;
        template<class VectorType, class Executor>
        void forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result) const noexcept;
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_short(const MDCellType& cell) const { return force<Executor>(cell); }
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_long(const MDCellType& cell) const { return Vector<ScalarType>(cell.getDOF(), 0); }
        void swap(VASPModel& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const WorkingDirType& getWorkingDir() const noexcept { return workingDir; }
        [[nodiscard]] size_t getNumParticle() const noexcept;
    private:
        ProcessFuture run_vasp() const;
        friend class Test;
    };

    template<class ScalarType>
    VASPModel<ScalarType>::VASPModel(std::string pathToVasp_,
                                     const char* pathToINCAR,
                                     const char* pathToPOTCAR,
                                     const char* pathToKPOINTS,
                                     Utils::Array<size_t> numOfEachType_,
                                     unsigned int numMPIProcess_)
            : pathToVasp(std::move(pathToVasp_))
            , workingDir("/tmp/tmpXXXXXX")
            , numOfEachType(std::move(numOfEachType_))
            , numMPIProcess(numMPIProcess_) {
        using namespace Physica::Utils;
        auto path = makePath("%s/log", workingDir.getName());
        logFd = open(path.data(), O_WRONLY | O_TRUNC | O_CREAT, S_IRUSR | S_IWUSR);
        if (logFd == -1)
            throw SyscallException();

        path = makePath("%s/INCAR", workingDir.getName());
        copyFile(pathToINCAR, path.data());
        path = makePath("%s/POTCAR", workingDir.getName());
        copyFile(pathToPOTCAR, path.data());
        path = makePath("%s/KPOINTS", workingDir.getName());
        copyFile(pathToKPOINTS, path.data());
    }

    template<class ScalarType>
    template<class Executor>
    Vector<ScalarType> VASPModel<ScalarType>::force(const MDCellType& cell) const {
        Vector<ScalarType> result(cell.getNumParticle());
        forceAsync<Vector<ScalarType>, Executor>(cell, result);
        return result;
    }

    template<class ScalarType>
    template<class VectorType, class Executor>
    void VASPModel<ScalarType>::forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result) const noexcept {
        try {
            /* Make POSCAR */ {
                const auto path = Utils::makePath("%s/POSCAR", workingDir.getName());
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
            const auto path = Utils::makePath("%s/OUTCAR", workingDir.getName());
            const Outcar outcar(path.data(), getNumParticle());
            result.getDerived() = outcar.getForce();
        }
        catch (std::exception& e) {
            fprintf(stderr, "%s\n", e.what());
            exit(EXIT_FAILURE);
        }
    }

    template<class ScalarType>
    void VASPModel<ScalarType>::swap(VASPModel& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        pathToVasp.swap(obj.pathToVasp);
        workingDir.swap(obj.workingDir);
        std::swap(logFd, obj.logFd);
        numOfEachType.swap(obj.numOfEachType);
        std::swap(numMPIProcess, obj.numMPIProcess);
    }

    template<class ScalarType>
    size_t VASPModel<ScalarType>::getNumParticle() const noexcept {
        size_t result = 0;
        for (size_t elem : numOfEachType)
            result += elem;
        return result;
    }

    template<class ScalarType>
    ProcessFuture VASPModel<ScalarType>::run_vasp() const {
        return ProcessExecutor::schedule([&]() {
            const int standardErr = dup(STDERR_FILENO);
            dup2(logFd, STDOUT_FILENO);
            dup2(logFd, STDERR_FILENO);
            [[maybe_unused]] int err = chdir(workingDir.getName());
            /* Execute */ {
                constexpr int bufferLength = 16; // 16 is enough for unsigned int
                char numProcess[bufferLength];
                [[maybe_unused]] int count = sprintf(numProcess, "%d", numMPIProcess);
                assert(0 <= count && count < bufferLength && "[Error]: Unexpected bad printf");
                execlp("mpirun", "mpirun", "-n", numProcess, pathToVasp.c_str(), nullptr);
            }
            dup2(standardErr, STDERR_FILENO);
            perror("[Error]: Failed to execute VASP");
            _exit(EXIT_FAILURE);
        });
    }
}
