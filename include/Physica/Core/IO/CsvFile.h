/*
 * Copyright 2023-2024 Weibo He.
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
#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h>

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
        using DataTypeArray = Array<DataType>;
        using DefaultValueArray = Array<std::optional<DefaultValue>>;
    private:
        using ArrayPointer = void*;

        DataTypeArray datatypes;
        Array<std::string> headers;
        Array<ArrayPointer> data;

        Array<std::optional<DefaultValue>> defaultValues;
    public:
        CsvFile(DataTypeArray datatypes_, const char* path);
        CsvFile(DataTypeArray datatypes_, Array<std::optional<DefaultValue>> defaultValues_, const char* path);
        CsvFile(const CsvFile&) = default;
        CsvFile(CsvFile&&) noexcept = default;
        ~CsvFile();
        /* Operators */
        CsvFile& operator=(CsvFile obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<Scalar T>
        [[nodiscard]] VectorND<T> asVector(size_t col) const;
        size_t countMissingValue(size_t col, const char* missingTag = "NA") const;
        void swap(CsvFile& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const DataTypeArray& getDatatypes() const noexcept { return datatypes; }
        [[nodiscard]] const Array<std::string>& getHeaders() const noexcept { return headers; }
        [[nodiscard]] inline const Array<signed char>& asChars(size_t col) const noexcept;
        [[nodiscard]] inline const Array<unsigned char>& asUsignedChars(size_t col) const noexcept;
        [[nodiscard]] inline const Array<short>& asShorts(size_t col) const noexcept;
        [[nodiscard]] inline const Array<unsigned short>& asUnsignedShorts(size_t col) const noexcept;
        [[nodiscard]] inline const Array<int>& asInts(size_t col) const noexcept;
        [[nodiscard]] inline const Array<unsigned int>& asUnsignedInts(size_t col) const noexcept;
        [[nodiscard]] inline const Array<long>& asLongs(size_t col) const noexcept;
        [[nodiscard]] inline const Array<unsigned long>& asUnsignedLongs(size_t col) const noexcept;
        [[nodiscard]] inline const VectorND<float32>& asFloats(size_t col) const noexcept;
        [[nodiscard]] inline const VectorND<float64>& asDoubles(size_t col) const noexcept;
        [[nodiscard]] inline const Array<bool>& asBools(size_t col) const noexcept;
        [[nodiscard]] inline const Array<std::string>& asStrings(size_t col) const noexcept;
        [[nodiscard]] inline size_t getRow() const noexcept;
        [[nodiscard]] size_t getCol() const noexcept { return datatypes.getLength(); }
        [[nodiscard]] const Array<std::optional<DefaultValue>>& getDefaultValues() const noexcept { return defaultValues; }
        /* Static members */
        [[nodiscard]] static const char* toString(DataType type);
    private:
        void allocate();
        void readFile(const char* path);
    };

    template<Scalar T>
    VectorND<T> CsvFile::asVector(size_t col) const {
        const DataType type = datatypes[col];
        VectorND<T> result(getRow());
        if (type == FLOAT)
            result = asFloats(col);
        else if (type == DOUBLE)
            result = asDoubles(col);
        else {
            for (size_t i = 0; i < result.getLength(); ++i) {
                switch (type) {
                case CHAR:
                    result[i] = T(asChars(col)[i]);
                    break;
                case UCHAR:
                    result[i] = T(asUsignedChars(col)[i]);
                    break;
                case SHORT:
                    result[i] = T(asShorts(col)[i]);
                    break;
                case USHORT:
                    result[i] = T(asUnsignedShorts(col)[i]);
                    break;
                case INT:
                    result[i] = T(asInts(col)[i]);
                    break;
                case UINT:
                    result[i] = T(asUnsignedInts(col)[i]);
                    break;
                case LONG:
                    result[i] = T(asLongs(col)[i]);
                    break;
                case ULONG:
                    result[i] = T(asUnsignedLongs(col)[i]);
                    break;
                case BOOL:
                    result[i] = T(asBools(col)[i]);
                    break;
                default: [[unlikely]]
                    throw std::invalid_argument("[Error]: Cannot convert string type to float type");
                }
            }
        }
        return result;
    }

    inline const Array<signed char>& CsvFile::asChars(size_t col) const noexcept {
        assert(datatypes[col] == CHAR && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Array<signed char>*>(data[col]);
    }

    inline const Array<unsigned char>& CsvFile::asUsignedChars(size_t col) const noexcept {
        assert(datatypes[col] == UCHAR && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Array<unsigned char>*>(data[col]);
    }

    inline const Array<short>& CsvFile::asShorts(size_t col) const noexcept {
        assert(datatypes[col] == SHORT && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Array<short>*>(data[col]);
    }

    inline const Array<unsigned short>& CsvFile::asUnsignedShorts(size_t col) const noexcept {
        assert(datatypes[col] == USHORT && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Array<unsigned short>*>(data[col]);
    }

    inline const Array<int>& CsvFile::asInts(size_t col) const noexcept {
        assert(datatypes[col] == INT && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Array<int>*>(data[col]);
    }

    inline const Array<unsigned int>& CsvFile::asUnsignedInts(size_t col) const noexcept {
        assert(datatypes[col] == UINT && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Array<unsigned int>*>(data[col]);
    }

    inline const Array<long>& CsvFile::asLongs(size_t col) const noexcept {
        assert(datatypes[col] == LONG && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Array<long>*>(data[col]);
    }

    inline const Array<unsigned long>& CsvFile::asUnsignedLongs(size_t col) const noexcept {
        assert(datatypes[col] == ULONG && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Array<unsigned long>*>(data[col]);
    }

    inline const VectorND<float32>& CsvFile::asFloats(size_t col) const noexcept {
        assert(datatypes[col] == FLOAT && "[Error]: Invalid conversion");
        return *reinterpret_cast<const VectorND<float32>*>(data[col]);
    }

    inline const VectorND<float64>& CsvFile::asDoubles(size_t col) const noexcept {
        assert(datatypes[col] == DOUBLE && "[Error]: Invalid conversion");
        return *reinterpret_cast<const VectorND<float64>*>(data[col]);
    }

    inline const Array<bool>& CsvFile::asBools(size_t col) const noexcept {
        assert(datatypes[col] == BOOL && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Array<bool>*>(data[col]);
    }

    inline const Array<std::string>& CsvFile::asStrings(size_t col) const noexcept {
        assert(datatypes[col] == STRING && "[Error]: Invalid conversion");
        return *reinterpret_cast<const Array<std::string>*>(data[col]);
    }

    inline size_t CsvFile::getRow() const noexcept {
        //Use char here but any datatype is ok if we need length only
        return reinterpret_cast<Array<char>*>(data[0])->getLength();
    }
}
