#ifndef ALIGNED_H
#define ALIGNED_H

#include <xmmintrin.h>

template <typename T, size_t O, size_t A>
class Aligned {
public:
    Aligned(size_t size) {
	    m_v = reinterpret_cast<T *>(_mm_malloc((size+O)*sizeof(T), A)) + O;
    }
    ~Aligned() { if (m_v) _mm_free(m_v - O); }
    // Conversion to pointer
    operator T*() { return m_v; }
    // Array indexing
    float &operator[](int i) { return m_v[i]; }
    const float &operator[](int i) const { return m_v[i]; }
    // Forbid normal copy
    Aligned(const Aligned<T,O,A> &other) = delete;
    // Forbid normal assignment
    Aligned<T,O,A> &operator=(Aligned<T,O,A> &other) = delete;
    // Move copy constructor operator
    Aligned(Aligned<T,O,A> &&other) {
        m_v = other.m_v;
        other.m_v = nullptr;
    }
    // Move assignment operator
    Aligned<T,O,A> &operator=(Aligned<T,O,A> &&other) {
        m_v = other.m_v;
        other.m_v = nullptr;
        return *this;
    }
private:
    float *m_v;
};

#endif
