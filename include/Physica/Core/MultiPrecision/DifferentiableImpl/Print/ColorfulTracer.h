/*
 * Copyright 2024 Weibo He.
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

#include "Physica/Utils/Unix/ColorGuard.h"

namespace Physica::Core {
    template<class ScalarType, unsigned int Order> class DiffTracer;

    template<class ScalarType, unsigned int Order>
    class ColorfulTracer {
        using TracerType = DiffTracer<ScalarType, Order>;

        const TracerType& tracer;
    public:
        ColorfulTracer(const TracerType& tracer_);
        ColorfulTracer(const ColorfulTracer&) = delete;
        ColorfulTracer(ColorfulTracer&&) noexcept = delete;
        ~ColorfulTracer() = default;
        /* Operators */
        ColorfulTracer& operator=(const ColorfulTracer&) = delete;
        ColorfulTracer& operator=(ColorfulTracer&&) noexcept = delete;
        template<class AnyScalar, unsigned int AnyOrder>
        friend std::ostream& operator<<(std::ostream& os, const ColorfulTracer<AnyScalar, AnyOrder>& obj);
    };

    template<class ScalarType, unsigned int Order>
    ColorfulTracer<ScalarType, Order>::ColorfulTracer(const TracerType& tracer_) : tracer(tracer_) {}

    template<class AnyScalar, unsigned int AnyOrder>
    std::ostream& operator<<(std::ostream& os, const ColorfulTracer<AnyScalar, AnyOrder>& obj) {
        using Color = typename Utils::ColorGuard::Color;
        using ColorGuard = Utils::ColorGuard;
        const auto& tracer = obj.tracer;
        const auto& list = tracer.getTraceList();
        for (auto ite = list.cbegin(); ite != list.cend(); ++ite) {
            const auto& tracer = *ite;
            {
                ColorGuard guard(os, Color::Yellow, false);
                os << "Trace " << tracer.getValues().data() << ":\n";
            }
            os << tracer.color() << '\n';
        }
        return os;
    }
}
