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
#include <fstream>
#include <iostream>
#include "Physica/Core/IO/CsvFile.h"
#include "Physica/Core/Exception/IOException.h"
#include "Physica/Core/Exception/BadFileFormatException.h"

namespace Physica::Core {
    CsvFile::CsvFile(DataTypeArray datatypes_, const char* path)
            : datatypes(std::move(datatypes_)) {
        allocate();
        readFile(path);
    }

    CsvFile::CsvFile(DataTypeArray datatypes_, Utils::Array<std::optional<DefaultValue>> defaultValues_, const char* path)
            : datatypes(std::move(datatypes_))
            , defaultValues(std::move(defaultValues_)) {
        assert(defaultValues.getLength() == getColumn() && "[Error]: Not enough values");
        allocate();
        readFile(path);
    }

    CsvFile::~CsvFile() {
        for (size_t i = 0; i < getColumn(); ++i) {
            using Physica::Utils::Array;
            void* p_array = data[i];
            switch (datatypes[i]) {
            case CHAR:
                delete reinterpret_cast<Array<char>*>(p_array);
                break;
            case UCHAR:
                delete reinterpret_cast<Array<unsigned char>*>(p_array);
                break;
            case SHORT:
                delete reinterpret_cast<Array<short>*>(p_array);
                break;
            case USHORT:
                delete reinterpret_cast<Array<unsigned short>*>(p_array);
                break;
            case INT:
                delete reinterpret_cast<Array<int>*>(p_array);
                break;
            case UINT:
                delete reinterpret_cast<Array<unsigned int>*>(p_array);
                break;
            case LONG:
                delete reinterpret_cast<Array<long>*>(p_array);
                break;
            case ULONG:
                delete reinterpret_cast<Array<unsigned long>*>(p_array);
                break;
            case FLOAT:
                delete reinterpret_cast<Vector<float32>*>(p_array);
                break;
            case DOUBLE:
                delete reinterpret_cast<Vector<float64>*>(p_array);
                break;
            case BOOL:
                delete reinterpret_cast<Array<bool>*>(p_array);
                break;
            case STRING:
                delete reinterpret_cast<Array<std::string>*>(p_array);
            default: [[unlikely]]
                assert("This is impossible");
            }
        }
    }

    size_t CsvFile::countMissingValue(size_t column, const char* missingTag) const {
        const auto datatype = datatypes[column];
        if (datatype != DataType::STRING && !defaultValues[column].has_value())
            return 0;

        size_t result = 0;
        switch (datatype) {
        case CHAR: {
            const signed char value = defaultValues[column]->char_value;
            const auto& arr = asChars(column);
            for (auto c : arr)
                result += value == c;
            break;
        }
        case UCHAR: {
            const unsigned char value = defaultValues[column]->uchar_value;
            const auto& arr = asUsignedChars(column);
            for (auto c : arr)
                result += value == c;
            break;
        }
        case SHORT: {
            const short value = defaultValues[column]->short_value;
            const auto& arr = asShorts(column);
            for (auto s : arr)
                result += value == s;
            break;
        }
        case USHORT: {
            const unsigned short value = defaultValues[column]->ushort_value;
            const auto& arr = asUnsignedShorts(column);
            for (auto s : arr)
                result += value == s;
            break;
        }
        case INT: {
            const int value = defaultValues[column]->int_value;
            const auto& arr = asInts(column);
            for (auto i : arr)
                result += value == i;
            break;
        }
        case UINT: {
            const unsigned int value = defaultValues[column]->uint_value;
            const auto& arr = asUnsignedInts(column);
            for (auto i : arr)
                result += value == i;
            break;
        }
        case LONG: {
            const long value = defaultValues[column]->long_value;
            const auto& arr = asLongs(column);
            for (auto l : arr)
                result += value == l;
            break;
        }
        case ULONG: {
            const unsigned long value = defaultValues[column]->ulong_value;
            const auto& arr = asUnsignedLongs(column);
            for (auto l : arr)
                result += value == l;
            break;
        }
        case FLOAT: {
            const float32 value = defaultValues[column]->float_value;
            const auto& arr = asFloats(column);
            for (auto f : arr)
                result += value == f;
            break;
        }
        case DOUBLE: {
            const float64 value = defaultValues[column]->double_value;
            const auto& arr = asDoubles(column);
            for (auto d : arr)
                result += value == d;
            break;
        }
        case BOOL: {
            const bool value = defaultValues[column]->bool_value;
            const auto& arr = asBools(column);
            for (auto b : arr)
                result += value == b;
            break;
        }
        case STRING: {
            const auto& arr = asStrings(column);
            for (const auto& s : arr)
                result += s == missingTag;
            break;
        }
        default: [[unlikely]]
            assert("[Error]: This is impossible");
        }
        return result;
    }

    void CsvFile::swap(CsvFile& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        datatypes.swap(obj.datatypes);
        headers.swap(obj.headers);
        data.swap(obj.data);
        defaultValues.swap(obj.defaultValues);
    }

    const char* CsvFile::toString(DataType type) {
        switch (type) {
        case CHAR:
            return "signed char";
        case UCHAR:
            return "unsigned char";
        case SHORT:
            return "short";
        case USHORT:
            return "unsigned short";
        case INT:
            return "int";
        case UINT:
            return "unsigned int";
        case LONG:
            return "long";
        case ULONG:
            return "unsigned long";
        case FLOAT:
            return "float";
        case DOUBLE:
            return "double";
        case BOOL:
            return "bool";
        case STRING:
            return "string";     
        default: [[unlikely]]
            throw std::invalid_argument("[Error]: Unrecognized DataType");
        }
    }

    void CsvFile::allocate() {
        headers.resize(getColumn());
        defaultValues.resize(getColumn());
        data.resize(getColumn());
        for (size_t i = 0; i < getColumn(); ++i) {
            using Physica::Utils::Array;
            void* p_array;
            switch (datatypes[i]) {
            case CHAR:
                p_array = new Array<signed char>();
                break;
            case UCHAR:
                p_array = new Array<unsigned char>();
                break;
            case SHORT:
                p_array = new Array<short>();
                break;
            case USHORT:
                p_array = new Array<unsigned short>();
                break;
            case INT:
                p_array = new Array<int>();
                break;
            case UINT:
                p_array = new Array<unsigned int>();
                break;
            case LONG:
                p_array = new Array<long>();
                break;
            case ULONG:
                p_array = new Array<unsigned long>();
                break;
            case FLOAT:
                p_array = new Vector<float32>();
                break;
            case DOUBLE:
                p_array = new Vector<float64>();
                break;
            case BOOL:
                p_array = new Array<bool>();
                break;
            case STRING:
                p_array = new Array<std::string>();
                break;        
            default: [[unlikely]]
                throw std::invalid_argument("[Error]: Unrecognized DataType");
            }
            data[i] = p_array;
        }
    }

    void CsvFile::readFile(const char* path) {
        std::ifstream fin(path);
        if (!fin)
            throw IOException("[Error]: File does not exist");
        std::string buffer{};
        for (size_t i = 0; i < getColumn(); ++i) {
            const bool isLast = i == getColumn() - 1;
            std::getline(fin, buffer, isLast ? '\n' : ',');
            headers[i] = buffer;
        }

        while (fin.good()) {
            {
                char c = fin.peek();
                while (std::isspace(c)) {
                    fin.get();
                    c = fin.peek();
                }
                if (c == EOF)
                    break;
            }

            for (size_t column = 0; column < getColumn(); ++column) {
                using Physica::Utils::Array;
                const bool isLast = column == getColumn() - 1;
                void* const p_array = data[column];
                const auto& defaultValue = defaultValues[column];
                switch (datatypes[column]) {
                case CHAR: {
                    int c; //Use int instread of char to avoid specialization for character
                    fin >> c;
                    if (!fin && defaultValue.has_value()) {
                        fin.clear();
                        std::getline(fin, buffer, isLast ? '\n' : ',');
                        c = defaultValue->char_value;
                    }
                    reinterpret_cast<Array<char>*>(p_array)->append(c);
                    break;
                }
                case UCHAR: {
                    unsigned int c; //Use unsigned int instread of unsigned char to avoid specialization for character
                    fin >> c;
                    if (!fin && defaultValue.has_value()) {
                        fin.clear();
                        std::getline(fin, buffer, isLast ? '\n' : ',');
                        c = defaultValue->char_value;
                    }
                    reinterpret_cast<Array<char>*>(p_array)->append(c);
                    break;
                }
                case SHORT: {
                    short s;
                    fin >> s;
                    if (!fin && defaultValue.has_value()) {
                        fin.clear();
                        std::getline(fin, buffer, isLast ? '\n' : ',');
                        s = defaultValue->short_value;
                    }
                    reinterpret_cast<Array<short>*>(p_array)->append(s);
                    break;
                }
                case USHORT: {
                    unsigned short s;
                    fin >> s;
                    if (!fin && defaultValue.has_value()) {
                        fin.clear();
                        std::getline(fin, buffer, isLast ? '\n' : ',');
                        s = defaultValue->ushort_value;
                    }
                    reinterpret_cast<Array<unsigned short>*>(p_array)->append(s);
                    break;
                }
                case INT: {
                    int i;
                    fin >> i;
                    if (!fin && defaultValue.has_value()) {
                        fin.clear();
                        std::getline(fin, buffer, isLast ? '\n' : ',');
                        i = defaultValue->int_value;
                    }
                    reinterpret_cast<Array<int>*>(p_array)->append(i);
                    break;
                }
                case UINT: {
                    unsigned int i;
                    fin >> i;
                    if (!fin && defaultValue.has_value()) {
                        fin.clear();
                        std::getline(fin, buffer, isLast ? '\n' : ',');
                        i = defaultValue->uint_value;
                    }
                    reinterpret_cast<Array<unsigned int>*>(p_array)->append(i);
                    break;
                }
                case LONG: {
                    long l;
                    fin >> l;
                    if (!fin && defaultValue.has_value()) {
                        fin.clear();
                        std::getline(fin, buffer, isLast ? '\n' : ',');
                        l = defaultValue->long_value;
                    }
                    reinterpret_cast<Array<long>*>(p_array)->append(l);
                    break;
                }
                case ULONG: {
                    unsigned long l;
                    fin >> l;
                    if (!fin && defaultValue.has_value()) {
                        fin.clear();
                        std::getline(fin, buffer, isLast ? '\n' : ',');
                        l = defaultValue->ulong_value;
                    }
                    reinterpret_cast<Array<unsigned long>*>(p_array)->append(l);
                    break;
                }
                case FLOAT: {
                    float32 f;
                    fin >> f;
                    if (!fin && defaultValue.has_value()) {
                        fin.clear();
                        std::getline(fin, buffer, isLast ? '\n' : ',');
                        f = defaultValue->float_value;
                    }
                    reinterpret_cast<Vector<float32>*>(p_array)->append(f);
                    break;
                }
                case DOUBLE: {
                    float64 d;
                    fin >> d;
                    if (!fin && defaultValue.has_value()) {
                        fin.clear();
                        std::getline(fin, buffer, isLast ? '\n' : ',');
                        d = defaultValue->double_value;
                    }
                    reinterpret_cast<Vector<float64>*>(p_array)->append(d);
                    break;
                }
                case BOOL: {
                    bool b;
                    fin >> b;
                    if (!fin && defaultValue.has_value()) {
                        fin.clear();
                        std::getline(fin, buffer, isLast ? '\n' : ',');
                        b = defaultValue->bool_value;
                    }
                    reinterpret_cast<Array<bool>*>(p_array)->append(b);
                    break;
                }
                case STRING: {
                    std::getline(fin, buffer, isLast ? '\n' : ',');
                    reinterpret_cast<Array<std::string>*>(p_array)->append(buffer);
                    break;
                }
                default: [[unlikely]]
                    throw std::invalid_argument("[Error]: Unrecognized DataType");
                }

                if (fin.peek() == ',')
                    fin.get();
                if (!fin)
                    throw BadFileFormatException("[Error]: Bad csv file");
            }
        }

        if (!fin)
            throw BadFileFormatException("[Error]: Bad csv file");
    }
}
