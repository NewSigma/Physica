/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Utils/Container/Array.h"
#include "Physica/Core/Parallel/Algorithm/Thread.h"
#include <stack>
#include "Test.h"

using namespace Physica;

namespace {
    template<class C>
    consteval void rangeTest() noexcept {
        using I = std::ranges::iterator_t<C>;
        static_assert(std::indirectly_readable<I>);
        static_assert(std::indirectly_writable<I, typename C::value_type>);
        static_assert(std::incrementable<I>);
        static_assert(std::sized_sentinel_for<I, I>);
        static_assert(std::contiguous_iterator<I>);

        static_assert(std::ranges::sized_range<C>);
        static_assert(std::ranges::contiguous_range<C>);
        static_assert(std::ranges::common_range<C>);
        static_assert(std::ranges::viewable_range<C>);
    }

    void structuredBinding() {
        Array<long, 3> arr{1, 2, 3};
        auto [x, y, z] = arr;
        expect(x == 1);
        expect(y == 2);
        expect(z == 3);
    }

    void emptyCopy() {
        // Test that we allow copying an empty array. This is useful when we default-initialize members of an object.
        Array<int> arr{};
        Array<int> copy = arr;
        expect(copy.getCapacity() == 0);
    }

    void makeStack() noexcept {
        // Test that std::stack accepts Array
        std::stack<int, Array<int>> s;
        s.push(1);
        expect(s.top() == 1 && s.size() == 1);
        s.pop();
        expect(s.empty());
    }

    void array_read() {
        Array<int, 3> origin{1, 2, 3};
        auto fixed = Array<int, 3>::read(origin.data());
        auto dynamic = Array<int>::read(3, origin.data());
        expect(fixed == origin);
        expect(dynamic == Array<int>(origin));
    }

    void generate() {
        constexpr size_t N = 32;
        /* Empty */ {
            auto arr = Array<int>::generate([](size_t) { return 0; }, 0);
            expect(arr.getLength() == 0);
        }
        /* Parallizable */ {
            auto fixed = Array<size_t, 32>::generate<Thread>([](size_t i) { return i; });
            auto dyn = Array<size_t>::generate<Thread>([](size_t i) { return i; }, N);
            expect(dyn.getLength() == N);
            for (size_t i = 0; i < N; ++i)
                expect(fixed[i] == i && dyn[i] == i);
        }
        /* Non-default-constructable */ {
            struct NoDefault {
                size_t value;

                explicit NoDefault(size_t value_) : value(value_) {}
            };
            auto arr = Array<NoDefault>::generate([](size_t i) { return NoDefault(i); }, N);
            expect(arr.getLength() == N);
            for (size_t i = 0; i < N; ++i)
                expect(arr[i].value == i);
        }
    }

    void compare() {
        Array<int> a, b;
        expect(a == b);
    }
}

int main() {
    rangeTest<Array<long, 3>>();
    rangeTest<Array<long>>();
    structuredBinding();
    emptyCopy();
    makeStack();
    array_read();
    generate();
    compare();
    return 0;
}
