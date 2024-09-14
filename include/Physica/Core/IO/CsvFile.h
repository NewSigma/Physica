/*
 * Copyright 2023 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

namespace Physica::Core {
    class PHYSICA_API CsvFile {
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
            float32 float_value;
            float64 double_value;
            bool bool_value;

            DefaultValue(signed char c) : char_value(c) {}
            DefaultValue(unsigned char c) : uchar_value(c) {}
            DefaultValue(short s) : short_value(s) {}
            DefaultValue(unsigned short s) : ushort_value(s) {}
            DefaultValue(int i) : int_value(i) {}
            DefaultValue(unsigned int i) : uint_value(i) {}
            DefaultValue(long l) : long_value(l) {}
            DefaultValue(unsigned long l) : ulong_value(l) {}
            DefaultValue(float32 f) : float_value(std::move(f)) {}
            DefaultValue(float64 d) : double_value(std::move(d)) {}
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
        template<class ScalarType>
        [[nodiscard]] Vector<ScalarType> asVector(size_t column) const;
        size_t countMissingValue(size_t column, const char* missingTag = "NA") const;
        void swap(CsvFile& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const DataTypeArray& getDatatypes() const noexcept { return datatypes; }
        [[nodiscard]] size_t getColumn() const noexcept { return datatypes.getLength(); }
        [[nodiscard]] const Utils::Array<std::string>& getHeaders() const noexcept { return headers; }
        [[nodiscard]] inline const Utils::Array<signed char>& asChars(size_t column) const noexcept;
        [[nodiscard]] inline const Utils::Array<unsigned char>& asUsignedChars(size_t column) const noexcept;
        [[nodiscard]] inline const Utils::Array<short>& asShorts(size_t column) const noexcept;
        [[nodiscard]] inline const Utils::Array<unsigned short>& asUnsignedShorts(size_t column) const noexcept;
        [[nodiscard]] inline const Utils::Array<int>& asInts(size_t column) const noexcept;
        [[nodiscard]] inline const Utils::Array<unsigned int>& asUnsignedInts(size_t column) const noexcept;
        [[nodiscard]] inline const Utils::Array<long>& asLongs(size_t column) const noexcept;
        [[nodiscard]] inline const Utils::Array<unsigned long>& asUnsignedLongs(size_t column) const noexcept;
        [[nodiscard]] inline const Vector<float32>& asFloats(size_t column) const noexcept;
        [[nodiscard]] inline const Vector<float64>& asDoubles(size_t column) const noexcept;
        [[nodiscard]] inline const Utils::Array<bool>& asBools(size_t column) const noexcept;
        [[nodiscard]] inline const Utils::Array<std::string>& asStrings(size_t column) const noexcept;
        [[nodiscard]] inline size_t getRow() const noexcept;
        [[nodiscard]] const Utils::Array<std::optional<DefaultValue>>& getDefaultValues() const noexcept { return defaultValues; }
        /* Static members */
        [[nodiscard]] static const char* toString(DataType type);
    private:
        void allocate();
        void readFile(const char* path);
    };

    template<class ScalarType>
    Vector<ScalarType> CsvFile::asVector(size_t column) const {
        const DataType type = datatypes[column];
        Vector<ScalarType> result(getRow());
        if (type == FLOAT)
            result = asFloats(column);
        else if (type == DOUBLE)
            result = asDoubles(column);
        else {
            for (size_t i = 0; i < result.getLength(); ++i) {
                switch (type) {
                case CHAR:
                    result[i] = ScalarType(asChars(column)[i]);
                    break;
                case UCHAR:
                    result[i] = ScalarType(asUsignedChars(column)[i]);
                    break;
                case SHORT:
                    result[i] = ScalarType(asShorts(column)[i]);
                    break;
                case USHORT:
                    result[i] = ScalarType(asUnsignedShorts(column)[i]);
                    break;
                case INT:
                    result[i] = ScalarType(asInts(column)[i]);
                    break;
                case UINT:
                    result[i] = ScalarType(asUnsignedInts(column)[i]);
                    break;
                case LONG:
                    result[i] = ScalarType(asLongs(column)[i]);
                    break;
                case ULONG:
                    result[i] = ScalarType(asUnsignedLongs(column)[i]);
                    break;
                case BOOL:
                    result[i] = ScalarType(asBools(column)[i]);
                    break;
                default: [[unlikely]]
                    throw std::invalid_argument("[Error]: Cannot convert string type to float type");
                }
            }
        }
        return result;
    }

    inline const Utils::Array<signed char>& CsvFile::asChars(size_t column) const noexcept {
        assert(datatypes[column] == CHAR && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Utils::Array<signed char>*>(data[column]);
    }

    inline const Utils::Array<unsigned char>& CsvFile::asUsignedChars(size_t column) const noexcept {
        assert(datatypes[column] == UCHAR && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Utils::Array<unsigned char>*>(data[column]);
    }

    inline const Utils::Array<short>& CsvFile::asShorts(size_t column) const noexcept {
        assert(datatypes[column] == SHORT && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Utils::Array<short>*>(data[column]);
    }

    inline const Utils::Array<unsigned short>& CsvFile::asUnsignedShorts(size_t column) const noexcept {
        assert(datatypes[column] == USHORT && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Utils::Array<unsigned short>*>(data[column]);
    }

    inline const Utils::Array<int>& CsvFile::asInts(size_t column) const noexcept {
        assert(datatypes[column] == INT && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Utils::Array<int>*>(data[column]);
    }

    inline const Utils::Array<unsigned int>& CsvFile::asUnsignedInts(size_t column) const noexcept {
        assert(datatypes[column] == UINT && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Utils::Array<unsigned int>*>(data[column]);
    }

    inline const Utils::Array<long>& CsvFile::asLongs(size_t column) const noexcept {
        assert(datatypes[column] == LONG && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Utils::Array<long>*>(data[column]);
    }

    inline const Utils::Array<unsigned long>& CsvFile::asUnsignedLongs(size_t column) const noexcept {
        assert(datatypes[column] == ULONG && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Utils::Array<unsigned long>*>(data[column]);
    }

    inline const Vector<float32>& CsvFile::asFloats(size_t column) const noexcept {
        assert(datatypes[column] == FLOAT && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Vector<float32>*>(data[column]);
    }

    inline const Vector<float64>& CsvFile::asDoubles(size_t column) const noexcept {
        assert(datatypes[column] == DOUBLE && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Vector<float64>*>(data[column]);
    }

    inline const Utils::Array<bool>& CsvFile::asBools(size_t column) const noexcept {
        assert(datatypes[column] == BOOL && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Utils::Array<bool>*>(data[column]);
    }

    inline const Utils::Array<std::string>& CsvFile::asStrings(size_t column) const noexcept {
        assert(datatypes[column] == STRING && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Utils::Array<std::string>*>(data[column]);
    }

    inline size_t CsvFile::getRow() const noexcept {
        //Use char here but any datatype is ok if we need length only
        return reinterpret_cast<Utils::Array<char>*>(data[0])->getLength();
    }
}
