# Public interface of Array

There are three specialization of template Array, they all share a set of interface:

1. operator[]
2. empty()
3. Iterator
4. operator== operator!=
5. init()
    Initialize a element at given position, double init() will lead to memory leak.(Except fixed array)
6. clear()
    This function cooperates with init(), when we have to remove all elements and reinit them all.
7. getLength()
8. setLength()
9. getCapacity()

Notice that 5, 6, 8 have different meanings between two specialization.
