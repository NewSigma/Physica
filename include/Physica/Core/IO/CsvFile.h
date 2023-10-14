/*
 * Copyright 2023 WeiBo He.
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

#include <optional>
#include <variant>
#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Utils/Container/Array/Array.h"

namespace Physica::Core {
    class CsvFile {
    public:
        enum DataType {
            CHAR,
            UCHAR,
            SHORT,
            USHORT,
            INT,
            UINT,
            LONG,
            ULONG,
            FLOAT,
            DOUBLE,
            BOOL,
            STRING
        };
 
        union DefaultValue {
            signed char char_value;
            unsigned char uchar_value;
            short short_value;
            unsigned short ushort_value;
            int int_value;
            unsigned int uint_value;
            long long_value;
            unsigned long ulong_value;
            float float_value;
            double double_value;
            bool bool_value;

            DefaultValue(signed char c) : char_value(c) {}
            DefaultValue(unsigned char c) : uchar_value(c) {}
            DefaultValue(short s) : short_value(s) {}
            DefaultValue(unsigned short s) : ushort_value(s) {}
            DefaultValue(int i) : int_value(i) {}
            DefaultValue(unsigned int i) : uint_value(i) {}
            DefaultValue(long l) : long_value(l) {}
            DefaultValue(unsigned long l) : ulong_value(l) {}
            DefaultValue(float f) : float_value(f) {}
            DefaultValue(double d) : double_value(d) {}
            DefaultValue(bool b) : bool_value(b) {}
        };
        using DataTypeArray = Utils::Array<DataType>;
        using DefaultValueArray = Utils::Array<std::optional<DefaultValue>>;
    private:
        using ArrayPointer = void*;

        DataTypeArray datatypes;
        Utils::Array<std::string> headers;
        Utils::Array<ArrayPointer> data;

        Utils::Array<std::optional<DefaultValue>> defaultValues;
    public:
        CsvFile(DataTypeArray datatypes_, const char* path);
        CsvFile(DataTypeArray datatypes_, Utils::Array<std::optional<DefaultValue>> defaultValues_, const char* path);
        CsvFile(const CsvFile&) = default;
        CsvFile(CsvFile&&) noexcept = default;
        ~CsvFile();
        /* Operators */
        CsvFile& operator=(CsvFile obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(CsvFile& obj) noexcept;
        /* Getters */
        [[nodiscard]] const DataTypeArray& getDatatypes() const noexcept { return datatypes; }
        [[nodiscard]] size_t getColumn() const noexcept { return datatypes.getLength(); }
        [[nodiscard]] inline size_t getRow() const noexcept;
        /* Static members */
        [[nodiscard]] static const char* toString(DataType type);
    private:
        void allocate();
        void readFile(const char* path);
    };

    inline size_t CsvFile::getRow() const noexcept {
        //Use char here but any datatype is ok if we need length only
        return reinterpret_cast<Utils::Array<char>*>(data[0])->getLength();
    }
}
