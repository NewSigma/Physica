/*
 * Copyright 2023-2026 Weibo He.
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
#include "Physica/Core/IO/HDF5/HDF5.h"
#include "Physica/Core/Utils/Unix/TempFile.h"
#include "Test.h"

using namespace Physica;

namespace {
    void expect_throw(auto&& fn) noexcept {
        bool threw = false;
        try {
            fn();
        }
        catch (...) {
            threw = true;
        }
        expect(threw);
    }

    void predicates() {
        TempFile temp("/tmp/tmpXXXXXX");
        auto h5f = H5File::open(temp.getName());
        auto group = H5Group::create(h5f, "G");
        const auto type = H5Type::get<int>();
        const auto space = H5DataSpace<1>(1);
        auto dataset = h5f.createDataSet<1>("D", type, space);
        auto attr = h5f.createAttribute("A", type, space);

        expect(h5f.isa<H5File>() && !h5f.isa<H5Group>());
        expect(group.isa<H5Group>());
        expect(type.isa<H5Type>());
        expect(space.isa<H5DataSpace<1>>());
        expect(dataset.isa<H5Dataset<1>>());
        expect(attr.isa<H5Attribute>());

        auto casted = std::move(group).cast<H5Group>();
        expect(casted.isValid() && casted.isa<H5Group>());
    }

    void stringIO() {
        TempFile temp("/tmp/tmpXXXXXX");
        const char* str = "This is a str\nstr line2";

        const auto dataspace = H5DataSpace<1>(strlen(str));
        {
            auto h5f = H5File::open(temp.getName());
            auto dataset = h5f.createDataSet<1>("/data", H5Type::get<char>(), dataspace);
            dataset.write(str, H5Type::get<char>());
        }
        {
            Array<char, 32> buffer{};
            auto h5f = H5File::open(temp.getName(), H5File::ReadOnly);
            auto dataset = h5f.openDataSet<1>("/data");
            dataset.readStr(buffer.data());
            expect(strcmp(str, buffer.data()) == 0);
            expect(h5f.isReadOnly());
            expect(dataset.isReadOnly());
        }
    }

    void readOnlyWrite() {
        TempFile temp("/tmp/tmpXXXXXX");
        const auto type = H5Type::get<int>();
        const auto space = H5DataSpace<1>(1);
        const int value = 42;
        {
            auto h5f = H5File::open(temp.getName());
            auto dataset = h5f.createDataSet<1>("/data", type, space);
            dataset.write(&value, type);
            auto attr = h5f.createAttribute("A", type, space);
            attr.write(type, &value);
        }

        auto h5f = H5File::open(temp.getName(), H5File::ReadOnly);

        expect_throw([&] { std::ignore = h5f.createDataSet<1>("/new", type, space); });
        expect_throw([&] { std::ignore = h5f.createAttribute("B", type, space); });
        expect_throw([&] { std::ignore = H5Group::create(h5f, "G"); });
        expect_throw([&] { h5f.openDataSet<1>("/data").write(&value, type); });
        expect_throw([&] { h5f.openAttribute("A").write(type, &value); });
    }
}

int main() {
    predicates();
    stringIO();
    readOnlyWrite();
    return 0;
}
