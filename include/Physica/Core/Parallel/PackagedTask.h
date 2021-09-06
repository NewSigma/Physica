/*
 * Copyright 2021 WeiBo He.
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

namespace Physica::Core::Parallel {
    class Task {
    public:
        virtual ~Task() = default;
        virtual void execute() = 0;
    };

    template<class Function, class... Args>
    class PackagedTask : public Task {
        std::packaged_task<Function(Args...)> task;
    public:
        PackagedTask(std::packaged_task<Function(Args...)> task_) : task(std::move(task_)) {}
        ~PackagedTask() override = default;
        /* Operations */
        void execute() override { task(); }
    };
}
