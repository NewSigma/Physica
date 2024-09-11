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
    template<class ScalarType, unsigned int Order> class TraceSegment;

    template<class ScalarType, unsigned int Order>
    class ColorfulTraceSegment {
        using SegmentType = TraceSegment<ScalarType, Order>;

        const SegmentType& segment;
    public:
        ColorfulTraceSegment(const SegmentType& segment_);
        ColorfulTraceSegment(const ColorfulTraceSegment&) = delete;
        ColorfulTraceSegment(ColorfulTraceSegment&&) noexcept = delete;
        ~ColorfulTraceSegment() = default;
        /* Operators */
        ColorfulTraceSegment& operator=(const ColorfulTraceSegment&) = delete;
        ColorfulTraceSegment& operator=(ColorfulTraceSegment&&) noexcept = delete;
        template<class AnyScalar, unsigned int AnyOrder>
        friend std::ostream& operator<<(std::ostream& os, const ColorfulTraceSegment<AnyScalar, AnyOrder>& obj);
    };

    template<class ScalarType, unsigned int Order>
    ColorfulTraceSegment<ScalarType, Order>::ColorfulTraceSegment(const SegmentType& segment_) : segment(segment_) {}

    template<class AnyScalar, unsigned int AnyOrder>
    std::ostream& operator<<(std::ostream& os, const ColorfulTraceSegment<AnyScalar, AnyOrder>& obj) {
        using Color = typename Utils::ColorGuard::Color;
        using ColorGuard = Utils::ColorGuard;
        const auto& segment = obj.segment;
        const size_t length = segment.getLength();
        for (size_t i = 0; i < length; ++i) {
            const auto source = segment.getRecords()[i].source;
            const auto width = os.width();
            /* Print ID */ {
                ColorGuard guard(os, Color::Green, true);
                os << std::setw(4) << i << std::setw(width) << ": ";
            }
            {
                ColorGuard guard(os, Color::Magenta, true);
                os << std::setw(10) << ExprTypeToStr(source) << std::setw(width) << ' ';
            }

            /* Print value */ {
                ColorGuard guard(os, Color::Cyan, true);
                os << segment.getValues()[i] << ' ' << segment.getGrads()[i] << ' ';
            }
            os << "Op: ";
            const size_t idFirstOperand = segment.getRecords()[i].idFirstOperand;
            const size_t num = TraceSegment<AnyScalar, AnyOrder>::numOperand(source);
            for (size_t j = 0; j < num; ++j) {
                const auto& op = segment.getOperands()[idFirstOperand + j];
                const size_t index = segment.find(op);
                const bool found = index < length;
                if (found) {
                    ColorGuard guard(os, Color::Green, true);
                    os << index;
                }
                else {
                    ColorGuard guard(os, Color::Yellow, false);
                    os << op.value_ptr();
                }
                os << ' ';
            }
            os << '\n';
        }
        return os;
    }
}
