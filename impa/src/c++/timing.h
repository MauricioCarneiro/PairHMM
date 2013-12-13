#ifndef TIMING_H
#define TIMING_H

#include <chrono>
#include <iostream>

class Timing {       // Report elapsed lifetime of object
public:
    Timing(): m_start(now()) {}
    ~Timing() {
        using us = std::chrono::microseconds;
        double t = std::chrono::duration_cast<us>(now()-m_start).count();
        std::cerr << t/1.0e6 << "\n";
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
