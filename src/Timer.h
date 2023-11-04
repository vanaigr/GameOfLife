#pragma once
#include <chrono>

template<class Units = std::chrono::microseconds>
class Timer {
private:
    std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
public:
    Timer() = default;

    uint32_t elapsedTime() const {
        auto elapsedTime = std::chrono::duration_cast<Units>(std::chrono::steady_clock::now() - startTime).count();
        return elapsedTime;
    }
};
