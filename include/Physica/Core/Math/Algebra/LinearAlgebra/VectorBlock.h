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
namespace Physica::Core {
    template<class T>
    class VectorBlock {
    public:
        using ScalarType = typename T::ScalarType;
        using VectorType = T;
    private:
        T& vec;
        size_t from;
        size_t to;
    public:
        VectorBlock(T& vec_, size_t from_, size_t to_);
        VectorBlock(T& vec_, size_t from_);
        VectorBlock(const VectorBlock&) = delete;
        VectorBlock(VectorBlock&&) noexcept = delete;
        ~VectorBlock() = default;
        /* Operators */
        VectorBlock& operator=(const VectorBlock&) = delete;
        VectorBlock& operator=(VectorBlock&&) noexcept = delete;
        ScalarType& operator[](size_t index) { assert((index + from) < to); return vec[index + from]; }
        const ScalarType& operator[](size_t index) const { assert((index + from) < to); return vec[index + from]; }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return to - from; }
    };

    template<class T>
    VectorBlock<T>::VectorBlock(T& vec_, size_t from_, size_t to_) : vec(vec_), from(from_), to(to_) {
        assert(from_ < to_);
        assert(to_ <= vec.getLength());
    }

    template<class T>
    VectorBlock<T>::VectorBlock(T& vec_, size_t from_) : VectorBlock(vec_, from_) {}
}
