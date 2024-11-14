#ifndef ARENA_H
#define ARENA_H

#include <iostream>

class Arena {
public:
    void* allocate(size_t sz) {
        void* ptr = ::operator new(sz); // Simple allocation
        return ptr;
    }

    void deallocate(void* p) {
        ::operator delete(p);
    }
};

void* operator new(size_t sz, Arena& a) {
    return a.allocate(sz);
}

template<class T>
void destroy(T* p, Arena& a) {
    if (p) {
        p->~T();        // Call destructor explicitly
        a.deallocate(p); // Return memory to the arena
    }
}

#endif