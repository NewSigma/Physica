/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

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
int main() {
    int a = -1;
    //Verify two's complement.
    int b = ~a + 1;
    if(b != -a)
        return 1;
    //Verify cast between two's complement and unsigned dose not change the bits.
    unsigned int c = static_cast<unsigned int>(a);
    if((a & c) != a)
        return 1;
    return 0;
}