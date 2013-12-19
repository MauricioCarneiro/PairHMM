#ifndef TIMING_H
#define TIMING_H

#include <chrono>
#include <iostream>

class Timing {
public:
    double elapsed(void) {
        using us = std::chrono::microseconds;
        double t = std::chrono::duration_cast<us>(now()-m_start).count();
        return t/1.0e6;
    }
    void mark(void) {
        m_start = now();
    }
private:
    using time_point =
        std::chrono::time_point<std::chrono::high_resolution_clock>;
    static time_point now() {
        return std::chrono::high_resolution_clock::now();
    }
    time_point m_start;
};

#endif
